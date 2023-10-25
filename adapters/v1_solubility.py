from adapters import register_adapter
from pydantic import BaseModel
from typing import Optional
from wrappers.registry import get_wrapper_registry
from wrappers.legacy_solubility.default import (
    LegacySolubilityInput,
    LegacySolubilityResult,
    LegacySolubilityOutput,
    LegacySolubilityResponse
)

class V1SolubilityInput(BaseModel):
    solvent: str #(str): solvent SMILES
    solute: str #(str): solute SMILES
    temp: float #(float): temperature for solubility prediction (K)
    ref_solvent: Optional[str] = None #(str, optional): reference solvent SMILES
    ref_solubility: Optional[float] = None #(float, optional): solute solubility in reference solvent as logS (log10(mol/L))
    ref_temp: Optional[float] = None #(float, optional): temperature for reference solubility (K)
    hsub298: Optional[float] = None #(float, optional): sublimation enthalpy of solute at 298 K (kcal/mol)
    cp_gas_298: Optional[float] = None #(float, optional): gas phase heat capacity of solute at 298 K (cal/K/mol)
    cp_solid_298: Optional[float] = None #(float, optional): solid phase heat capacity of solute at 298 K (cal/K/mol)

class V1SolubilityLegacyInput(BaseModel):
    task_list: list[V1SolubilityInput]

class V1SolubilityAsyncReturn(BaseModel):
    request: V1SolubilityInput
    task_id: str


class V1SolubilityResult(BaseModel):
    __root__: list[LegacySolubilityResult]


@register_adapter(
    name="v1_solubility",
    input_class=V1SolubilityInput,
    result_class=V1SolubilityResult
)
class V1SolubilityAdapter:
    """V1 Solubility Adapter"""
    prefix = "legacy/solubility"
    methods = ["POST"]

    def call_sync(self, input: V1SolubilityInput) -> V1SolubilityResult:
        wrapper = get_wrapper_registry().get_wrapper(module="legacy_solubility")
        wrapper_input = self.convert_input(input=input)
        print(wrapper_input)
        wrapper_response = wrapper.call_sync(wrapper_input)
        print("Here is the out put")
        print(wrapper_response)
        response = self.convert_response(wrapper_response=wrapper_response)

        return response

    async def __call__(self, input: V1SolubilityInput, priority: int = 0
                       ) -> V1SolubilityAsyncReturn:
        from askcos2_celery.tasks import legacy_task

        async_result = legacy_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        async_return = V1SolubilityAsyncReturn(
            request=input, task_id=task_id)

        return async_return

    # def __call__(self, input: V1SolubilityInput) -> V1SolubilityResult:
    #     return self.call_sync(input)

    @staticmethod
    def convert_input(input: V1SolubilityInput) -> LegacySolubilityInput:
        wrapper_input = LegacySolubilityInput(task_list=[input])
        return wrapper_input

    @staticmethod
    def convert_response(wrapper_response: LegacySolubilityResponse
                         ) -> V1SolubilityResult:
        result = V1SolubilityResult(__root__=wrapper_response)

        return result
