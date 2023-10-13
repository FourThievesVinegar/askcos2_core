from pydantic import BaseModel
from schemas.base import LowerCamelAliasModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class AtomMapWLNInput(LowerCamelAliasModel):
    smiles: list[str]


class AtomMapWLNOutput(BaseModel):
    error: str
    status: str
    results: list[str]


class AtomMapWLNResponse(BaseResponse):
    result: list[str]


@register_wrapper(
    name="atom_map_wln",
    input_class=AtomMapWLNInput,
    output_class=AtomMapWLNOutput,
    response_class=AtomMapWLNResponse
)
class AtomMapWLNWrapper(BaseWrapper):
    """Wrapper class for WLN Atom Mapper"""
    prefixes = ["atom_map/wln"]

    def call_sync(self, input: AtomMapWLNInput) -> AtomMapWLNResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: AtomMapWLNInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapWLNResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: AtomMapWLNOutput
                                   ) -> AtomMapWLNResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = AtomMapWLNResponse(**response)

        return response
