from datetime import datetime
from fastapi import Depends
from pydantic import BaseModel, constr, Field
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
    cluster_method: str = "rxn_class"
    fp_type: str = "morgan"
    fp_length: int = 512
    fp_radius: int = 1
    classification_threshold: float = 0.2


class ExpandOneOptions(BaseModel):
    # aliasing to v1 fields
    template_max_count: int = Field(default=100, alias="template_count")
    template_max_cum_prob: int = Field(default=0.995, alias="max_cum_template_prob")
    banned_chemicals: list[str] = Field(default=[], alias="forbidden_molecules")
    banned_reactions: list[str] = Field(default=[], alias="known_bad_reactions")

    retro_backend_options: list[RetroBackendOption] = [RetroBackendOption()]
    use_fast_filter: bool = True
    filter_threshold: float = 0.75
    retro_rerank_backend: str | None = None
    cluster_precursors: bool = False
    cluster_setting: ClusterSetting = None
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
    max_trees: int | None = None
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
    json_format: Literal["treedata", "nodelink"] = "treedata"
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
    expand_one_options: ExpandOneOptions = ExpandOneOptions()
    build_tree_options: BuildTreeOptions = BuildTreeOptions()
    enumerate_paths_options: EnumeratePathsOptions = EnumeratePathsOptions()


class MCTSResult(BaseModel):
    stats: dict[str, Any] | None
    paths: list[dict[str, Any]]
    graph: dict[str, Any] | None
    version: int | None


class MCTSOutput(BaseModel):
    error: str
    status: str
    results: MCTSResult


class MCTSResponse(BaseResponse):
    result: MCTSResult


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

        # actual backend call
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

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
        from askcos2_celery.tasks import task_with_token
        async_result = task_with_token.apply_async(
            args=(self.name, input.dict(), token), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> MCTSResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: MCTSOutput
                                   ) -> MCTSResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = MCTSResponse(**response)

        return response
