from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel


class AtomMapIndigoInput(BaseModel):
    smiles: list[str]


class AtomMapIndigoOutput(BaseModel):
    error: str
    status: str
    results: list[list[str]]


class AtomMapIndigoResponse(BaseResponse):
    result: list[str]


@register_wrapper(
    name="atom_map_indigo",
    input_class=AtomMapIndigoInput,
    output_class=AtomMapIndigoOutput,
    response_class=AtomMapIndigoResponse
)
class AtomMapIndigoWrapper(BaseWrapper):
    """Wrapper class for Indigo Atom Mapper"""
    prefixes = ["atom_map/indigo"]

    def call_sync(self, input: AtomMapIndigoInput) -> AtomMapIndigoResponse:
        output = super().call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: AtomMapIndigoInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapIndigoResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: AtomMapIndigoOutput
                                   ) -> AtomMapIndigoResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result_format": "json",
            "result": output.results[0]
        }
        response = AtomMapIndigoResponse(**response)

        return response
