from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel


class FastFilterInput(BaseModel):
    smiles: list[str]


class FastFilterOutcome(BaseModel):
    smiles: str
    template_ids: list
    num_examples: int


class FastFilterResult(BaseModel):
    rank: float
    outcome: FastFilterOutcome
    score: float
    prob: float


class FastFilterOutput(BaseModel):
    error: str
    status: str
    results: list[list[list[FastFilterResult]]]


class FastFilterResponse(BaseResponse):
    result: FastFilterResult


@register_wrapper(
    name="fast_filter",
    input_class=FastFilterInput,
    output_class=FastFilterOutput,
    response_class=FastFilterResponse
)
class FastFilterWrapper(BaseWrapper):
    """Wrapper class for Fast Filter"""
    prefixes = ["fast_filter"]

    def call_sync(self, input: FastFilterInput) -> FastFilterResponse:
        if not self.prediction_url.endswith("fast_filter_evaluate"):
            self.prediction_url = f"{self.prediction_url}/fast_filter_evaluate"
            # this creates an infinite stacking loop w/o if statement
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: FastFilterInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> FastFilterResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: FastFilterOutput
                                   ) -> FastFilterResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results[0][0][0]
        }
        response = FastFilterResponse(**response)

        return response
