from wrappers import register_wrapper
from wrappers.base import BaseWrapper, BaseResponse
from pydantic import BaseModel
from typing import List, Dict


class ImpurityPredictorInput(BaseModel):
    rct_smi: str
    prd_smi: str
    sol_smi: str
    rea_smi: str


class ImpurityPredictorResult(BaseModel):
    predict_expand: list[Dict]


class ImpurityPredictorOutput(BaseModel):
    error: str
    status: str
    results: list[ImpurityPredictorResult]


class ImpurityPredictorResponse(BaseResponse):
    result: List[ImpurityPredictorResult]

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
        output = super().call_raw(input=input)
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

