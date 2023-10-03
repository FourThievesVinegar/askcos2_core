from pydantic import BaseModel
from wrappers import register_wrapper
from wrappers.base import BaseWrapper, BaseResponse


class PathwayRankerInput(BaseModel):
    trees: list[dict]
    clustering: bool = False
    cluster_method: str = "hdbscan"
    min_samples: int = 5
    min_cluster_size: int = 5


class PathwayRankerResult(BaseModel):
    scores: list
    encoded_trees: list
    clusters: list


class PathwayRankerOutput(BaseModel):
    error: str
    status: str
    results: list[PathwayRankerResult]


class PathwayRankerResponse(BaseResponse):
    result: PathwayRankerResult


@register_wrapper(
    name="pathway_ranker",
    input_class=PathwayRankerInput,
    output_class=PathwayRankerOutput,
    response_class=PathwayRankerResponse
)
class PathwayRankerWrapper(BaseWrapper):
    """Wrapper class for tree-LSTM based pathway ranker"""
    prefixes = ["pathway_ranker"]

    def call_sync(self, input: PathwayRankerInput) -> PathwayRankerResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: PathwayRankerInput, priority: int = 0) -> str:
        from askcos2_celery.tasks import pathway_ranker_task
        async_result = pathway_ranker_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> PathwayRankerResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: PathwayRankerOutput
                                   ) -> PathwayRankerResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results[0]
        }
        response = PathwayRankerResponse(**response)

        return response
