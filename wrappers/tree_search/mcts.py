import re
import traceback
import uuid
from datetime import datetime
from fastapi import Depends, HTTPException
from pydantic import BaseModel, Field, validator
from schemas.cluster import ClusterSetting
from schemas.retro import RetroBackendOption
from typing import Annotated, Any, Literal
from utils.oauth2 import oauth2_scheme
from utils.registry import get_util_registry
from utils.tree_search_results import TreeSearchSavedResults
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


def to_lower_camel(name: str) -> str:
    """
    Converts a snake_case string to lowerCamelCase
    """
    upper = "".join(word.capitalize() for word in name.split("_"))
    return upper[:1].lower() + upper[1:]


class ExpandOneOptions(BaseModel):
    # aliasing to v1 fields
    template_max_count: int = Field(default=100, alias="template_count")
    template_max_cum_prob: int = Field(default=0.995, alias="max_cum_template_prob")
    banned_chemicals: list[str] = Field(default_factory=list, alias="forbidden_molecules")
    banned_reactions: list[str] = Field(default_factory=list, alias="known_bad_reactions")

    retro_backend_options: list[RetroBackendOption] = [RetroBackendOption()]
    use_fast_filter: bool = True
    filter_threshold: float = 0.75
    retro_rerank_backend: Literal["relevance_heuristic", "scscore"] = None
    cluster_precursors: bool = False
    cluster_setting: ClusterSetting = ClusterSetting()
    extract_template: bool = False
    return_reacting_atoms: bool = True
    selectivity_check: bool = False


class BuildTreeOptions(BaseModel):
    expansion_time: int = 30
    max_iterations: int | None = None
    max_chemicals: int | None = None
    max_reactions: int | None = None
    max_templates: int | None = None
    max_branching: int = 25
    max_depth: int = 5
    exploration_weight: float = 1.0
    return_first: bool = False
    max_trees: int = 500
    max_ppg: float | None = None
    max_scscore: float | None = None
    max_elements: dict[str, int] | None = None
    min_history: dict[str, int] | None = None
    property_criteria: list[dict[str, Any]] | None = None
    termination_logic: dict[str, list[str]] = {"and": ["buyable"]}
    buyables_source: str | list[str] | None = None
    custom_buyables: list[str] | None = None


class EnumeratePathsOptions(BaseModel):
    path_format: Literal["json", "graph"] = "json"
    json_format: Literal["treedata", "nodelink"] = "nodelink"
    sorting_metric: Literal[
        "plausibility",
        "number_of_starting_materials",
        "number_of_reactions",
        "score"
    ] = "plausibility"
    validate_paths: bool = True
    score_trees: bool = False
    cluster_trees: bool = False
    cluster_method: Literal["hdbscan", "kmeans"] = "hdbscan"
    min_samples: int = 5
    min_cluster_size: int = 5
    paths_only: bool = False
    max_paths: int = 200


class MCTSInput(BaseModel):
    smiles: str
    description: str | None = ""
    tags: str | None = ""
    expand_one_options: ExpandOneOptions = ExpandOneOptions()
    build_tree_options: BuildTreeOptions = BuildTreeOptions()
    enumerate_paths_options: EnumeratePathsOptions = EnumeratePathsOptions()
    run_async: bool = False
    result_id: str = str(uuid.uuid4())

    @validator(*)


class MCTSResult(BaseModel):
    stats: dict[str, Any] | None
    paths: list[dict[str, Any]] | None
    graph: dict[str, Any] | None
    version: int | None = 2
    result_id: str = ""


class MCTSOutput(BaseModel):
    error: str
    status: str
    results: MCTSResult


class MCTSResponse(BaseResponse):
    result: MCTSResult


@register_wrapper(
    name="tree_search_mcts",
    input_class=MCTSInput,
    output_class=MCTSOutput,
    response_class=MCTSResponse
)
class MCTSWrapper(BaseWrapper):
    """Wrapper class for Monte Carlo Tree Search"""
    prefixes = ["tree_search/mcts"]
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
        input: MCTSInput,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> MCTSResponse:
        # banned_chemicals handling, requires login token
        if not input.expand_one_options.banned_chemicals:
            input.expand_one_options.banned_chemicals = []
        banned_chemicals_controller = get_util_registry().get_util(
            module="banned_chemicals"
        )
        user_banned_chemicals = banned_chemicals_controller.get(token=token).__root__
        user_banned_chemicals = [entry.smiles for entry in user_banned_chemicals
                                 if entry.active]
        input.expand_one_options.banned_chemicals.extend(user_banned_chemicals)

        # banned_chemicals handling, requires login token
        if not input.expand_one_options.banned_reactions:
            input.expand_one_options.banned_reactions = []
        banned_reactions_controller = get_util_registry().get_util(
            module="banned_reactions"
        )
        user_banned_reactions = banned_reactions_controller.get(token=token).__root__
        user_banned_reactions = [entry.smiles for entry in user_banned_reactions
                                 if entry.active]
        input.expand_one_options.banned_reactions.extend(user_banned_reactions)

        results_controller = get_util_registry().get_util(
            module="tree_search_results_controller"
        )

        if input.run_async:
            results_controller.update_result_state(
                result_id=input.result_id,
                state="started",
                token=token
            )
        try:
            # actual backend call
            output = self.call_raw(input=input)
            response = self.convert_output_to_response(output)
            result_doc = response.result
            result_doc.result_id = input.result_id
        except Exception:
            if input.run_async:
                results_controller.update_result_state(
                    result_id=input.result_id,
                    state="failed",
                    token=token
                )
            raise HTTPException(
                status_code=500,
                detail=f"mcts.call_sync() fails with the error: "
                       f"{traceback.format_exc()}"
            )

        if input.run_async:
            results_controller.update_result_state(
                result_id=input.result_id,
                state="completed",
                token=token
            )
            results_controller.save_results(
                result_id=input.result_id,
                result=result_doc,
                token=token
            )

        return response

    def call_sync_without_token(self, input: MCTSInput) -> MCTSResponse:
        """Special method to be called by other modules"""
        # banned_chemicals handling, does not require login token
        if not input.expand_one_options.banned_chemicals:
            input.expand_one_options.banned_chemicals = []

        # banned_chemicals handling, does not require login token
        if not input.expand_one_options.banned_reactions:
            input.expand_one_options.banned_reactions = []

        # actual backend call
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(
        self,
        input: MCTSInput,
        token: Annotated[str, Depends(oauth2_scheme)],
        priority: int = 0
    ) -> str:
        user_controller = get_util_registry().get_util(module="user_controller")
        user = user_controller.get_current_user(token)

        input.run_async = True
        input.result_id = str(uuid.uuid4())
        # Note that we can't use task_id as the result_id,
        # as it needs to be known beforehand
        settings = input.dict()
        settings = {k: v for k, v in settings.items() if "option" in k}

        saved_results = TreeSearchSavedResults(
            user=user.username,
            description=input.description,
            created=datetime.now(),
            modified=datetime.now(),
            result_id=input.result_id,
            result_state="pending",
            result_type="tree_builder",
            result=None,
            settings=settings,
            tags=[tag.strip() for tag in re.split(r"[,;]", input.tags)],
            shared_with=[user.username]
        )
        results_controller = get_util_registry().get_util(
            module="tree_search_results_controller"
        )
        results_controller.create(result=saved_results, token=token)

        from askcos2_celery.tasks import tree_search_mcts_task
        async_result = tree_search_mcts_task.apply_async(
            args=(self.name, input.dict(), token), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> MCTSResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: MCTSOutput
                                   ) -> MCTSResponse:
        status_code = 200
        message = "mcts.call_raw() successfully executed."
        result = output.results

        if not output.status == "SUCCESS":
            status_code = 500
            message = f"Backend error encountered during mcts.call_raw() " \
                      f"with the following error message {output.error}"

        response = MCTSResponse(
            status_code=status_code,
            message=message,
            result=result
        )

        return response
