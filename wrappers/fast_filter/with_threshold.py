from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel


class FastFilterWThresholdInput(BaseModel):
    smiles: list[str]
    threshold: float


class FastFilterWThresholdResult(BaseModel):
    flag: bool
    score: float


class FastFilterWThresholdOutput(BaseModel):
    error: str
    status: str
    results: list[FastFilterWThresholdResult]


class FastFilterWThresholdResponse(BaseResponse):
    result: FastFilterWThresholdResult


@register_wrapper(
    name="fast_filter_with_threshold",
    input_class=FastFilterWThresholdInput,
    output_class=FastFilterWThresholdOutput,
    response_class=FastFilterWThresholdResponse
)
class FastFilterWThresholdWrapper(BaseWrapper):
    """Wrapper class for Fast Filter with Threshold"""
    prefixes = ["fast_filter/with_threshold"]

    def call_sync(self, input: FastFilterWThresholdInput
                  ) -> FastFilterWThresholdResponse:
        if not self.prediction_url.endswith("filter_with_threshold"):
            self.prediction_url = f"{self.prediction_url}/filter_with_threshold"
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: FastFilterWThresholdInput, priority: int = 0
                         ) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> FastFilterWThresholdResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: FastFilterWThresholdOutput
                                   ) -> FastFilterWThresholdResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results[0]
        }
        response = FastFilterWThresholdResponse(**response)

        return response
