from pydantic import BaseModel
from schemas.base import LowerCamelAliasModel
from typing import Literal
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class ClusterInput(LowerCamelAliasModel):
    original: str
    outcomes: list[str]
    feature: Literal["original", "outcomes", "all"] = "original"
    cluster_method: Literal["kmeans", "hdbscan", "rxn_class"] = "kmeans"
    fp_type: Literal["morgan"] = "morgan"
    fp_length: int = 512
    fp_radius: int = 1
    scores: list[float] = None
    classification_threshold: float = 0.2


class ClusterOutput(BaseModel):
    error: str
    status: str
    results: list[list[int] | dict[str, str]]


class ClusterResponse(BaseResponse):
    result: list[list[int] | dict[str, str]]


@register_wrapper(
    name="cluster",
    input_class=ClusterInput,
    output_class=ClusterOutput,
    response_class=ClusterResponse
)
class ClusterWrapper(BaseWrapper):
    """Wrapper class for clustering reactions"""
    prefixes = ["cluster"]

    def call_sync(self, input: ClusterInput
                  ) -> ClusterResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: ClusterInput, priority: int = 0
                         ) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ClusterResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ClusterOutput
                                   ) -> ClusterResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = ClusterResponse(**response)

        return response
