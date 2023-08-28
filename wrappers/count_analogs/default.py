from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel
from typing import List

class CountAnalogsInput(BaseModel):
    smiles: list[str]


class CountAnalogsResult(BaseModel):
    __root__: int


class CountAnalogsOutput(BaseModel):
    error: str
    status: str
    results: int


class CountAnalogsResponse(BaseResponse):
    result: int


@register_wrapper(
    name="count_analogs",
    input_class=CountAnalogsInput,
    output_class=CountAnalogsOutput,
    response_class=CountAnalogsResponse
)
class CountAnalogsWrapper(BaseWrapper):
    """Wrapper class for Descriptors"""
    prefixes = ["count_analogs"]

    def call_sync(self, input: CountAnalogsInput) -> CountAnalogsResponse:
        output = super().call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: CountAnalogsInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> CountAnalogsResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: CountAnalogsOutput
                                   ) -> CountAnalogsResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = CountAnalogsResponse(**response)

        return response
