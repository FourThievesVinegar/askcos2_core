from pydantic import BaseModel
from wrappers import register_wrapper
from wrappers.base import BaseWrapper, BaseResponse


class PmiCalculatorInput(BaseModel):
    trees: list[dict]


class PmiCalculatorOutput(BaseModel):
    error: str
    status: str
    results: list[list[float]]


class PmiCalculatorResponse(BaseResponse):
    result: list[float]


@register_wrapper(
    name="pmi_calculator",
    input_class=PmiCalculatorInput,
    output_class=PmiCalculatorOutput,
    response_class=PmiCalculatorResponse
)
class PmiCalculatorWrapper(BaseWrapper):
    """Wrapper class for Pmi Calculator"""
    prefixes = ["pmi_calculator"]

    def call_sync(self, input: PmiCalculatorInput) -> PmiCalculatorResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)
        return response

    async def call_async(self, input: PmiCalculatorInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> PmiCalculatorResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: PmiCalculatorOutput
                                   ) -> PmiCalculatorResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results[0]
        }
        response = PmiCalculatorResponse(**response)

        return response
