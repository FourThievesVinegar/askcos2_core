from adapters import register_adapter
from pydantic import BaseModel
from typing import Optional, Literal
from wrappers.registry import get_wrapper_registry
from wrappers.forward.wldn5 import (
    ForwardWLDN5Input,
    ForwardWLDN5Result,
    ForwardWLDN5Response
)


class V1ForwardInput(BaseModel):
    reactants: str  #(str): SMILES string of reactants
    reagents: Optional[str] = None #(str, optional): SMILES string of reagents (default='')
    solvent: Optional[str] = None #(str, optional): SMILES string of solvent (default='')
    training_set: Optional[str] = "uspto_500k" #(str, optional): model training set to use
    model_version: Optional[Literal["wln_forward", "graph2smiles_forward_uspto_stereo",""]] = "" #(str, optional): model version to use for prediction (wln_forward (default) or graph2smiles_forward_uspto_stereo)
    num_results: Optional[int] = 100 #(int, optional): max number of results to return (default=100)
    atommap: Optional[bool] = False #(bool, optional): Flag to keep atom mapping from the prediction (default=False)
    priority: Optional[int] = 1 #(int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)


class V1ForwardAsyncReturn(BaseModel):
    request: V1ForwardInput
    task_id: str

class V1ForwardOutput(BaseModel):
    rank: int
    score: float
    prob: float
    mol_wt: float
    smiles: str


class V1ForwardResult(BaseModel):
    __root__: list[V1ForwardOutput]

@register_adapter(
    name="v1_forward",
    input_class=V1ForwardInput,
    result_class=V1ForwardResult
)
class V1ForwardAdapter:
    """V1 Forward Adapter"""
    prefix = "legacy/forward"
    methods = ["POST"]

    def call_sync(self, input: V1ForwardInput) -> V1ForwardResult:
        wrapper = get_wrapper_registry().get_wrapper(module="forward_wldn5")
        wrapper_input = self.convert_input(input=input)
        wrapper_response = wrapper.call_sync(wrapper_input)
        response = self.convert_response(wrapper_response=wrapper_response)

        return response

    async def __call__(self, input: V1ForwardInput, priority: int = 0
                       ) -> V1ForwardAsyncReturn:
        from askcos2_celery.tasks import legacy_task

        async_result = legacy_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        async_return = V1ForwardAsyncReturn(
            request=input, task_id=task_id)

        return async_return

    # def __call__(self, input: V1ForwardInput) -> V1ForwardResult:
    #     return self.call_sync(input)

    @staticmethod
    def convert_input(input: V1ForwardInput) -> ForwardWLDN5Input:
        wrapper_input = ForwardWLDN5Input(
            model_name=input.training_set,
            reactants=input.reactants,
            contexts=None
        )

        return wrapper_input

    @staticmethod
    def convert_response(wrapper_response: ForwardWLDN5Response
                         ) -> V1ForwardResult:
        results = wrapper_response.result[0]
        converted_results = [V1ForwardOutput(
            rank=result.rank,
            score=result.score,
            prob=result.prob,
            mol_wt=result.mol_wt,
            smiles=result.outcome.smiles
        )for result in results]
        final_result = V1ForwardResult(__root__=converted_results)

        return final_result
