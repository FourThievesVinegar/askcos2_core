from pydantic import BaseModel
from schemas.base import LowerCamelAliasModel
from typing import Literal
from wrappers import register_wrapper
from wrappers.atom_map.indigo import AtomMapIndigoInput, AtomMapIndigoResponse
from wrappers.atom_map.rxnmapper import AtomMapRXNMapperInput, AtomMapRXNMapperResponse
from wrappers.atom_map.wln import AtomMapWLNInput, AtomMapWLNResponse
from wrappers.base import BaseResponse, BaseWrapper
from wrappers.registry import get_wrapper_registry


class AtomMapInput(LowerCamelAliasModel):
    backend: Literal["indigo", "rxnmapper", "wln"] = "rxnmapper"
    smiles: list[str]


class AtomMapOutput(BaseModel):
    placeholder: str


class AtomMapResponse(BaseResponse):
    result: list[str]


@register_wrapper(
    name="atom_map_controller",
    input_class=AtomMapInput,
    output_class=AtomMapOutput,
    response_class=AtomMapResponse
)
class AtomMapController(BaseWrapper):
    """Atom Map Controller"""
    prefixes = ["atom_map/controller", "atom_map"]
    backend_wrapper_names = {
        "indigo": "atom_map_indigo",
        "rxnmapper": "atom_map_rxnmapper",
        "wln": "atom_map_wln"
    }

    def __init__(self):
        pass        # TODO: proper inheritance

    def call_sync(self, input: AtomMapInput) -> AtomMapResponse:
        module = self.backend_wrapper_names[input.backend]
        wrapper = get_wrapper_registry().get_wrapper(module=module)

        wrapper_input = self.convert_input(
            input=input, backend=input.backend)
        wrapper_response = wrapper.call_sync(wrapper_input)
        response = self.convert_response(
            wrapper_response=wrapper_response, backend=input.backend)

        return response

    async def call_async(self, input: AtomMapInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_input(
        input: AtomMapInput, backend: str
    ) -> AtomMapIndigoInput | AtomMapRXNMapperInput | AtomMapWLNInput:
        if backend == "indigo":
            wrapper_input = AtomMapIndigoInput(smiles=input.smiles)
        elif backend == "rxnmapper":
            wrapper_input = AtomMapRXNMapperInput(smiles=input.smiles)
        elif backend == "wln":
            wrapper_input = AtomMapWLNInput(smiles=input.smiles)
        else:
            raise ValueError(f"Unsupported atom map backend: {backend}!")

        return wrapper_input

    @staticmethod
    def convert_response(
        wrapper_response:
            AtomMapIndigoResponse | AtomMapRXNMapperResponse | AtomMapWLNResponse,
        backend: str
    ) -> AtomMapResponse:
        status_code = wrapper_response.status_code
        message = wrapper_response.message
        if backend in ["indigo", "wln"]:
            # list[str] -> list[str]
            result = wrapper_response.result
        elif backend == "rxnmapper":
            # list[RXNMapperResult] -> list[str]
            result = [r.mapped_rxn for r in wrapper_response.result]
        else:
            raise ValueError(f"Unsupported atom map backend: {backend}!")

        response = {
            "status_code": status_code,
            "message": message,
            "result": result
        }
        response = AtomMapResponse(**response)

        return response
