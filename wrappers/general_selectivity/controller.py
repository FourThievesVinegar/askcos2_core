from pydantic import BaseModel
from schemas.base import LowerCamelAliasModel
from typing import Literal
from wrappers import register_wrapper
from wrappers.general_selectivity.gnn import (
    GeneralSelectivityInput as GnnGeneralSelectivityInput,
    GeneralSelectivityResponse as GnnGeneralSelectivityResponse
)
from wrappers.general_selectivity.qm_gnn import (
    GeneralSelectivityInput as QmGnnGeneralSelectivityInput,
    GeneralSelectivityResponse as QmGnnGeneralSelectivityResponse
)
from wrappers.general_selectivity.qm_gnn_no_reagent import (
    GeneralSelectivityInput as NoReagentGeneralSelectivityInput,
    GeneralSelectivityResponse as NoReagentGeneralSelectivityResponse
)
from wrappers.base import BaseResponse, BaseWrapper
from wrappers.registry import get_wrapper_registry


class GeneralSelectivityInput(LowerCamelAliasModel):
    backend: Literal["gnn", "qm_gnn", "qm_gnn_no_reagent"] = "gnn"
    smiles: str
    atom_map_backend: Literal["indigo", "rxnmapper", "wln"] = "wln"
    mapped: bool = False
    all_outcomes: bool = False
    no_map_reagents: bool = False


class GeneralSelectivityOutput(BaseModel):
    placeholder: str


class GeneralSelectivityResult(BaseModel):
    smiles: str
    prob: float
    rank: int


class GeneralSelectivityResponse(BaseResponse):
    result: list[GeneralSelectivityResult]


@register_wrapper(
    name="general_selectivity_controller",
    input_class=GeneralSelectivityInput,
    output_class=GeneralSelectivityOutput,
    response_class=GeneralSelectivityResponse
)
class GeneralSelectivityController(BaseWrapper):
    """General Selectivity Controller"""
    prefixes = ["general_selectivity/controller", "general_selectivity"]
    backend_wrapper_names = {
        "gnn": "general_selectivity_gnn",
        "qm_gnn": "general_selectivity_qm_gnn",
        "qm_gnn_no_reagent": "general_selectivity_qm_gnn_no_reagent"
    }

    def __init__(self):
        pass        # TODO: proper inheritance

    def call_sync(self, input: GeneralSelectivityInput) -> GeneralSelectivityResponse:
        module = self.backend_wrapper_names[input.backend]
        wrapper = get_wrapper_registry().get_wrapper(module=module)

        wrapper_input = self.convert_input(
            input=input, backend=input.backend)
        wrapper_response = wrapper.call_sync(wrapper_input)
        response = self.convert_response(
            wrapper_response=wrapper_response, backend=input.backend)

        return response

    async def call_async(self, input: GeneralSelectivityInput, priority: int = 0
                         ) -> str:
        from askcos2_celery.tasks import general_selectivity_task
        async_result = general_selectivity_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> GeneralSelectivityResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_input(
        input: GeneralSelectivityInput, backend: str
    ) -> (GnnGeneralSelectivityInput | QmGnnGeneralSelectivityInput |
          NoReagentGeneralSelectivityInput):
        if backend == "gnn":
            wrapper_input = GnnGeneralSelectivityInput(
                smiles=input.smiles,
                atom_map_backend=input.atom_map_backend,
                mapped=input.mapped,
                all_outcomes=input.all_outcomes,
                no_map_reagents=input.no_map_reagents
            )
        elif backend == "qm_gnn":
            wrapper_input = QmGnnGeneralSelectivityInput(
                smiles=input.smiles,
                atom_map_backend=input.atom_map_backend,
                mapped=input.mapped,
                all_outcomes=input.all_outcomes,
                no_map_reagents=input.no_map_reagents
            )
        elif backend == "qm_gnn_no_reagent":
            wrapper_input = NoReagentGeneralSelectivityInput(
                smiles=input.smiles,
                atom_map_backend=input.atom_map_backend,
                mapped=input.mapped,
                all_outcomes=input.all_outcomes,
                no_map_reagents=input.no_map_reagents
            )
        else:
            raise ValueError(f"Unsupported atom map backend: {backend}!")

        return wrapper_input

    @staticmethod
    def convert_response(
        wrapper_response:
            GnnGeneralSelectivityResponse | QmGnnGeneralSelectivityResponse |
            NoReagentGeneralSelectivityResponse,
        backend: str
    ) -> GnnGeneralSelectivityResponse:
        status_code = wrapper_response.status_code
        message = wrapper_response.message
        result = wrapper_response.result

        response = {
            "status_code": status_code,
            "message": message,
            "result": result
        }
        response = GeneralSelectivityResponse(**response)

        return response
