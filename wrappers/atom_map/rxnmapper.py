from pydantic import BaseModel
from schemas.base import LowerCamelAliasModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class RXNMapperResult(BaseModel):
    mapped_rxn: str
    confidence: float


class AtomMapRXNMapperInput(LowerCamelAliasModel):
    smiles: list[str]


class AtomMapRXNMapperOutput(BaseModel):
    error: str
    status: str
    results: list[list[RXNMapperResult]]


class AtomMapRXNMapperResponse(BaseResponse):
    result: list[RXNMapperResult]


@register_wrapper(
    name="atom_map_rxnmapper",
    input_class=AtomMapRXNMapperInput,
    output_class=AtomMapRXNMapperOutput,
    response_class=AtomMapRXNMapperResponse
)
class AtomMapRXNMapperWrapper(BaseWrapper):
    """Wrapper class for IBM RXNMapper"""
    prefixes = ["atom_map/rxnmapper"]

    def call_sync(self, input: AtomMapRXNMapperInput) -> AtomMapRXNMapperResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: AtomMapRXNMapperInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapRXNMapperResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: AtomMapRXNMapperOutput
                                   ) -> AtomMapRXNMapperResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results[0]
        }
        response = AtomMapRXNMapperResponse(**response)

        return response
