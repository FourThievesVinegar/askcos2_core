from datetime import datetime
from fastapi import Depends
from pydantic import BaseModel, constr
from typing import Annotated, Any, Literal
from utils.oauth2 import oauth2_scheme
from utils.registry import get_util_registry
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class RetroBackendOption(BaseModel):
    retro_backend: str = "template_relevance"
    retro_model_name: str = "reaxys"
    max_num_templates: int = 100
    max_cum_prob: float = 0.995
    attribute_filter: list[dict[str, Any]] = []


class ClusterSetting(BaseModel):
    feature: str = "original"
    cluster_method: Literal["rxn_class", "hdbscan", "kmeans"] = "rxn_class"
    fp_type: str = "morgan"
    fp_length: int = 512
    fp_radius: int = 1
    classification_threshold: float = 0.2


class ExpandOneInput(BaseModel):
    smiles: str
    retro_backend_options: list[RetroBackendOption] = [RetroBackendOption()]
    banned_chemicals: list[str] = None
    banned_reactions: list[str] = None
    use_fast_filter: bool = True
    fast_filter_threshold: float = 0.75
    retro_rerank_backend: str = "relevance_heuristic"
    cluster_precursors: bool = True
    cluster_setting: ClusterSetting = ClusterSetting()
    extract_template: bool = False
    return_reacting_atoms: bool = True
    selectivity_check: bool = False


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
    result: list[RetroResult]


class BlacklistedEntry(BaseModel):
    id: str = ""
    user: str
    description: constr(max_length=1000) = None
    created: datetime = datetime.now()
    dt: constr(max_length=200) = None
    smiles: constr(max_length=5000)
    active: bool = True


class BlacklistedChemicals(BaseModel):
    __root__: list[BlacklistedEntry]


class BlacklistedReactions(BaseModel):
    __root__: list[BlacklistedEntry]


@register_wrapper(
    name="tree_search_expand_one",
    input_class=ExpandOneInput,
    output_class=ExpandOneOutput,
    response_class=ExpandOneResponse
)
class ExpandOneWrapper(BaseWrapper):
    """Wrapper class for clustering reactions"""
    prefixes = ["tree_search/expand_one"]

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
        response = self.convert_output_to_response(output)

        return response

    async def call_async(
        self,
        input: ExpandOneInput,
        token: Annotated[str, Depends(oauth2_scheme)],
        priority: int = 0
    ) -> str:
        from askcos2_celery.tasks import task_with_token
        async_result = task_with_token.apply_async(
            args=(self.name, input.dict(), token), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> ExpandOneResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ExpandOneOutput
                                   ) -> ExpandOneResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = ExpandOneResponse(**response)

        return response
