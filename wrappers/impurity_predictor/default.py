from pydantic import BaseModel
from typing import Literal
from wrappers import register_wrapper
from wrappers.base import BaseWrapper, BaseResponse


class ImpurityPredictorInput(BaseModel):
    predictor_backend: Literal[
        "augmented_transformer", "graph2smiles", "wldn5"] = "wldn5"
    predictor_model_name: str = "pistachio"
    inspector: str = "Reaxys inspector"
    atom_map_backend: Literal["indigo", "rxnmapper", "wln"]
    topn_outcome: int = 3
    insp_threshold: float = 0.2
    check_mapping: bool = True

    rct_smi: str
    prd_smi: str = ""
    sol_smi: str = ""
    rea_smi: str = ""


class ForwardResult(BaseModel):
    outcome: str
    score: float


class ImpurityResult(BaseModel):
    no: int
    prd_smiles: str
    prd_mw: float
    rct_rea_sol: list[dict]
    insp_scores: list[float]
    avg_insp_score: float
    similarity_to_major: float
    modes: list[int]


class ImpurityPredictorResult(BaseModel):
    predict_expand: list[ImpurityResult]
    predict_normal: list[ForwardResult]


class ImpurityPredictorOutput(BaseModel):
    error: str
    status: str
    results: ImpurityPredictorResult


class ImpurityPredictorResponse(BaseResponse):
    result: ImpurityPredictorResult


@register_wrapper(
    name="impurity_predictor",
    input_class=ImpurityPredictorInput,
    output_class=ImpurityPredictorOutput,
    response_class=ImpurityPredictorResponse
)
class ImpurityPredictorWrapper(BaseWrapper):
    """Wrapper class for Impurity Predictor"""
    prefixes = ["impurity_predictor"]

    def call_sync(self, input: ImpurityPredictorInput
                  ) -> ImpurityPredictorResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: ImpurityPredictorInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ImpurityPredictorResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ImpurityPredictorOutput
                                   ) -> ImpurityPredictorResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = ImpurityPredictorResponse(**response)

        return response

