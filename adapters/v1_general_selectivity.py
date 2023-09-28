from adapters import register_adapter
from pydantic import BaseModel
from typing import Optional, Literal
from wrappers.registry import get_wrapper_registry
from wrappers.general_selectivity.controller import (
    GeneralSelectivityInput,
    GeneralSelectivityResult,
    GeneralSelectivityResponse
)


class V1GeneralSelectivityInput(BaseModel):
    reactants: str
    product: str
    reagents: Optional[str] = None
    solvent: Optional[str] = None
    mapped: bool = False
    all_outcomes: bool = False
    verbose: bool = True
    mapper: Optional[Literal["Transformer", "WLN atom mapper"]] #Literal["indigo", "rxnmapper", "wln"] = "wln"
    no_map_reagents: bool = True
    mode: Optional[Literal["GNN", "qm_GNN"]] #Literal["gnn", "qm_gnn", "qm_gnn_no_reagent"] = "gnn"

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
        if input.mapper=="Transformer" or input.mapper==None:
            if input.mode=="GNN" or input.mode==None:
                wrapper_input = GeneralSelectivityInput(
                    backend="gnn",
                    smiles=input.reactants+">"+input.reagents+">"+input.product,
                    atom_map_backend="rxnmapper",
                    mapped=input.mapped,
                    all_outcomes=input.all_outcomes,
                    no_map_reagents=input.no_map_reagents
                )
            elif input.mode=="qm_GNN":
                wrapper_input = GeneralSelectivityInput(
                    backend="qm_gnn",
                    smiles=input.reactants+">"+input.reagents+">"+input.product,
                    atom_map_backend="rxnmapper",
                    mapped=input.mapped,
                    all_outcomes=input.all_outcomes,
                    no_map_reagents=input.no_map_reagents
                )
        elif input.mapper=="WLN atom mapper":
            if input.mode=="GNN":
                wrapper_input = GeneralSelectivityInput(
                    backend="gnn",
                    smiles=input.reactants+">"+input.reagents+">"+input.product,
                    atom_map_backend="wln",
                    mapped=input.mapped,
                    all_outcomes=input.all_outcomes,
                    no_map_reagents=input.no_map_reagents
                )
            elif input.mode=="qm_GNN":
                wrapper_input = GeneralSelectivityInput(
                    backend="qm_gnn",
                    smiles=input.reactants+">"+input.reagents+">"+input.product,
                    atom_map_backend="wln",
                    mapped=input.mapped,
                    all_outcomes=input.all_outcomes,
                    no_map_reagents=input.no_map_reagents
                )
        else:
            print("there is something wrong")
            wrapper_input = None
        return wrapper_input

    @staticmethod
    def convert_response(wrapper_response: GeneralSelectivityResponse
                         ) -> V1GeneralSelectivityResult:
        result = wrapper_response.result
        result = V1GeneralSelectivityResult(__root__=result)

        return result
