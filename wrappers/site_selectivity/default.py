from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel


class SiteSelectivityInput(BaseModel):
    smiles: list[str]


class SiteSelectivityResult(BaseModel):
    smiles: str
    index: int
    task: str
    atom_scores: list[float]


class SiteSelectivityOutput(BaseModel):
    error: str
    status: str
    results: list[list[SiteSelectivityResult]]


class SiteSelectivityResponse(BaseResponse):
    result: list[list[SiteSelectivityResult]]


@register_wrapper(
    name="site_selectivity",
    input_class=SiteSelectivityInput,
    output_class=SiteSelectivityOutput,
    response_class=SiteSelectivityResponse
)
class SiteSelectivityWrapper(BaseWrapper):
    """Wrapper class for Site Selectivity"""
    prefixes = ["site_selectivity"]

    def call_sync(self, input: SiteSelectivityInput) -> SiteSelectivityResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: SiteSelectivityInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> SiteSelectivityResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: SiteSelectivityOutput
                                   ) -> SiteSelectivityResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = SiteSelectivityResponse(**response)

        return response
