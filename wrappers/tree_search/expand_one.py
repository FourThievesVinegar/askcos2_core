from fastapi import Depends
from pydantic import BaseModel, Field
from schemas.base import LowerCamelAliasModel
from schemas.cluster import ClusterSetting
from schemas.retro import RetroBackendOption
from typing import Annotated, Any, Literal
from utils.oauth2 import oauth2_scheme
from utils.registry import get_util_registry
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class ExpandOneInput(LowerCamelAliasModel):
    smiles: str
    retro_backend_options: list[RetroBackendOption] = [RetroBackendOption()]
    banned_chemicals: list[str] = Field(default_factory=list)
    banned_reactions: list[str] = Field(default_factory=list)
    use_fast_filter: bool = True
    fast_filter_threshold: float = 0.75
    retro_rerank_backend: Literal["relevance_heuristic", "scscore"] | None = "scscore"
    cluster_precursors: bool = False
    cluster_setting: ClusterSetting = ClusterSetting()
    extract_template: bool = False
    return_reacting_atoms: bool = True
    selectivity_check: bool = False
    group_by_strategy: bool = False


class RetroResult(BaseModel):
    # from retro_controller
    outcome: str
    model_score: float
    normalized_model_score: float
    template: dict[str, Any] | None

    # extended from postprocessing in expand_one_controller
    retro_backend: str
    retro_model_name: str
    plausibility: float | None
    rms_molwt: float
    num_rings: int
    scscore: float
    group_id: int | None
    group_name: str | None
    mapped_smiles: str | None
    reacting_atoms: list[int] | None
    selec_error: bool | None
    mapped_outcomes: str | None
    mapped_precursors: str | None

    score: float
    rank: int


class ExpandOneOutput(BaseModel):
    error: str
    status: str
    results: list[RetroResult]


class ExpandOneResponse(BaseResponse):
    result: list[RetroResult] | dict[int, list[RetroResult]]


@register_wrapper(
    name="tree_search_expand_one",
    input_class=ExpandOneInput,
    output_class=ExpandOneOutput,
    response_class=ExpandOneResponse
)
class ExpandOneWrapper(BaseWrapper):
    """Wrapper class for ones tep expansion"""
    prefixes = ["tree_search/expand_one"]
    methods_to_bind: dict[str, list[str]] = {
        "get_config": ["GET"],
        "get_doc": ["GET"],
        "call_sync": ["POST"],
        "call_sync_without_token": ["POST"],
        "call_async": ["POST"],
        "retrieve": ["GET"]
    }

    def call_sync(
        self,
        input: ExpandOneInput,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> ExpandOneResponse:
        # banned_chemicals handling, requires login token
        if not input.banned_chemicals:
            input.banned_chemicals = []
        banned_chemicals_controller = get_util_registry().get_util(
            module="banned_chemicals"
        )
        user_banned_chemicals = banned_chemicals_controller.get(token=token).__root__
        user_banned_chemicals = [entry.smiles for entry in user_banned_chemicals
                                 if entry.active]
        input.banned_chemicals.extend(user_banned_chemicals)

        # banned_chemicals handling, requires login token
        if not input.banned_reactions:
            input.banned_reactions = []
        banned_reactions_controller = get_util_registry().get_util(
            module="banned_reactions"
        )
        user_banned_reactions = banned_reactions_controller.get(token=token).__root__
        user_banned_reactions = [entry.smiles for entry in user_banned_reactions
                                 if entry.active]
        input.banned_reactions.extend(user_banned_reactions)

        # actual backend call
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output, input)

        return response

    def call_sync_without_token(self, input: ExpandOneInput) -> ExpandOneResponse:
        """Special method to be called by other tree searchers"""
        # banned_chemicals handling, does not require login token
        if not input.banned_chemicals:
            input.banned_chemicals = []

        # banned_chemicals handling, does not require login token
        if not input.banned_reactions:
            input.banned_reactions = []

        # actual backend call
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output, input)

        return response

    async def call_async(
        self,
        input: ExpandOneInput,
        token: Annotated[str, Depends(oauth2_scheme)],
        priority: int = 0
    ) -> str:
        from askcos2_celery.tasks import tree_search_expand_one_task
        async_result = tree_search_expand_one_task.apply_async(
            args=(self.name, input.dict(), token), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> ExpandOneResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(
        output: ExpandOneOutput,
        input: ExpandOneInput
    ) -> ExpandOneResponse:
        status_code = 200
        message = "expand_one.call_raw() successfully executed."
        result = output.results

        if not output.status == "SUCCESS":
            status_code = 500
            message = f"Backend error encountered during expand_one.call_raw() " \
                      f"with the following error message {output.error}"

        if input.group_by_strategy:
            grouped_results = {}
            for i, strategy in enumerate(input.retro_backend_options):
                grouped_result = []
                for prediction in result:
                    if (
                        prediction.retro_backend == strategy.retro_backend and
                        prediction.retro_model_name == strategy.retro_model_name
                    ):
                        grouped_result.append(prediction)
                grouped_results[i] = grouped_result
            result = grouped_results

        response = ExpandOneResponse(
            status_code=status_code,
            message=message,
            result=result
        )

        return response
