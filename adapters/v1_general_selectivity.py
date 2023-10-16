from adapters import register_adapter
from fastapi import HTTPException
from pydantic import BaseModel
from typing import Literal
from wrappers.registry import get_wrapper_registry
from wrappers.general_selectivity.controller import (
    GeneralSelectivityInput,
    GeneralSelectivityResult,
    GeneralSelectivityResponse
)


class V1GeneralSelectivityInput(BaseModel):
    reactants: str
    product: str
    reagents: str = ""
    solvent: str = ""
    mapped: bool = False
    all_outcomes: bool = False
    verbose: bool = True
    mapper: Literal["Transformer", "WLN atom mapper"] = "Transformer"
    no_map_reagents: bool = True
    mode: Literal["GNN", "qm_GNN"] = "GNN"


class V1GeneralSelectivityAsyncReturn(BaseModel):
    request: V1GeneralSelectivityInput
    task_id: str


class V1GeneralSelectivityResult(BaseModel):
    __root__: list[GeneralSelectivityResult]


@register_adapter(
    name="v1_general_selectivity",
    input_class=V1GeneralSelectivityInput,
    result_class=V1GeneralSelectivityResult
)
class V1GeneralSelectivityAdapter:
    """V1 General Selectivity Adapter"""
    prefix = "legacy/general_selectivity"
    methods = ["POST"]

    def call_sync(self, input: V1GeneralSelectivityInput) -> V1GeneralSelectivityResult:
        wrapper = get_wrapper_registry().get_wrapper(module="general_selectivity_controller")
        print("I am inside the call Sync now!!!!!!!!!!!!!")
        print(input)
        wrapper_input = self.convert_input(input=input)
        print("Here is the wrapper input:::::::::::::::::")
        print(wrapper_input)
        wrapper_response = wrapper.call_sync(wrapper_input)
        response = self.convert_response(wrapper_response=wrapper_response)

        return response

    async def __call__(self, input: V1GeneralSelectivityInput, priority: int = 0
                       ) -> V1GeneralSelectivityAsyncReturn:
        from askcos2_celery.tasks import legacy_task

        async_result = legacy_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        async_return = V1GeneralSelectivityAsyncReturn(
            request=input, task_id=task_id)

        return async_return

    @staticmethod
    def convert_input(input: V1GeneralSelectivityInput) -> GeneralSelectivityInput:
        if input.mapper == "Transformer" or input.mapper is None:
            atom_map_backend = "rxnmapper"
        elif input.mapper == "WLN atom mapper":
            atom_map_backend = "wln"
        else:
            raise ValueError(f"Unsupported atom_map_backend {input.mapper}")

        if input.mode == "GNN" or input.mode is None:
            backend = "gnn"
        elif input.mode == "qm_GNN":
            backend = "qm_gnn"
        else:
            raise ValueError(f"Unsupported mode {input.mode}")

        wrapper_input = GeneralSelectivityInput(
            backend=backend,
            smiles=input.reactants + ">" + input.reagents + ">" + input.product,
            atom_map_backend=atom_map_backend,
            mapped=input.mapped,
            all_outcomes=input.all_outcomes,
            no_map_reagents=input.no_map_reagents
        )

        return wrapper_input

    @staticmethod
    def convert_response(wrapper_response: GeneralSelectivityResponse
                         ) -> V1GeneralSelectivityResult:
        if not wrapper_response.status_code == 200:
            # reverse parse the exception in string format
            e = eval(wrapper_response.message.split("\n")[-3].strip())
            raise ValueError(e)

        result = wrapper_response.result
        result = V1GeneralSelectivityResult(__root__=result)

        return result
