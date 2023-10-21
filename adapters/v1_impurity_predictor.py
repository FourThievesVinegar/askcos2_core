from adapters import register_adapter
from pydantic import BaseModel
from typing import Optional
from wrappers.registry import get_wrapper_registry
from wrappers.impurity_predictor.default import (
    ImpurityPredictorInput,
    ImpurityPredictorResult,
    ImpurityPredictorResponse,
    ImpurityResult
)


class V1ImpurityPredictorInput(BaseModel):
    reactants: str #(str): SMILES string of reactants
    reagents: Optional[str] = "" #(str, optional): SMILES string of reagents
    products: Optional[str] = "" #(str, optional): SMILES string of products
    solvent: Optional[str] = "" #(str, optional): SMILES string of solvent
    top_k: Optional[int] = 3  #(int, optional): max number of results to return
    threshold: Optional[float] = 0.1 #(float, optional): probability threshold
    predictor: Optional[str] = "WLN forward predictor" #(str, optional): forward predictor to use
    inspector: Optional[str] = "WLN forward inspector" #(str, optional): reaction scorer to use
    mapper: Optional[str] = "WLN atom mapper" #(str, optional): reaction atom mapper to use
    training_set: Optional[str] = "uspto_500k" #(str, optional): WLN forward model training set to use
    model_version: Optional[str] = "" #(str, optional): WLN forward model version to use for prediction
    check_mapping: Optional[bool] = True #(bool, optional): whether to check atom mapping


class V1ImpurityPredictorAsyncReturn(BaseModel):
    request: V1ImpurityPredictorInput
    task_id: str

class ForwardOutcome(BaseModel):
    smiles: str

class ForwardResult(BaseModel):
    rank: int
    outcome: ForwardOutcome
    score: float


class V1ImpurityPredictorResult(BaseModel):
    predict_expand: list[ImpurityResult]
    predict_normal: list[list[ForwardResult]]


@register_adapter(
    name="v1_impurity_predictor",
    input_class=V1ImpurityPredictorInput,
    result_class=V1ImpurityPredictorResult
)
class V1ImpurityPredictorAdapter:
    """V1 ImpurityPredictor Adapter"""
    prefix = "legacy/impurity"
    methods = ["POST"]

    def call_sync(self, input: V1ImpurityPredictorInput) -> V1ImpurityPredictorResult:
        wrapper = get_wrapper_registry().get_wrapper(module="impurity_predictor")
        wrapper_input = self.convert_input(input=input)
        wrapper_response = wrapper.call_sync(wrapper_input)
        response = self.convert_response(wrapper_response=wrapper_response)

        return response

    async def __call__(self, input: V1ImpurityPredictorInput, priority: int = 0
                       ) -> V1ImpurityPredictorAsyncReturn:
        from askcos2_celery.tasks import legacy_task

        async_result = legacy_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        async_return = V1ImpurityPredictorAsyncReturn(
            request=input, task_id=task_id)

        return async_return


    @staticmethod
    def convert_input(input: V1ImpurityPredictorInput) -> ImpurityPredictorInput:
        wrapper_input = ImpurityPredictorInput(
            rct_smi=input.reactants,
            predictor_model_name=input.training_set,
            insp_threshold=input.threshold,
            inspector=input.inspector,
            atom_map_backend="wln"
        )

        return wrapper_input

    @staticmethod
    def convert_response(wrapper_response: ImpurityPredictorResponse
                         ) -> V1ImpurityPredictorResult:
        result = wrapper_response.result
        converted_results = []
        for idx, result in enumerate(result.predict_normal):
            forward_result = ForwardResult(
                rank=idx+1,
                outcome=ForwardOutcome(smiles=result.outcome),
                score=result.score
            )
            converted_results.append(forward_result)
        result = wrapper_response.result
        result = V1ImpurityPredictorResult(
            predict_expand = result.predict_expand,
            predict_normal = [converted_results]
        )

        return result
