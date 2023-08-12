from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from typing import List, Dict, Optional, Literal
from pydantic import BaseModel


class EvaluateReactionsInput(BaseModel):
    reactant_smiles: str
    target: str
    contexts: List[List]
    forward_scorer: Literal["Fast_Filter", "Template_Free"] = "Template_Free"
    top_n: Optional[int] = 10
    return_all_outcomes: Optional[bool]=False
    isomeric_smiles: Optional[bool]=False


class EvaluateReactionsResult(BaseModel):
    __root__: Dict


class EvaluateReactionsOutput(BaseModel):
    error: str
    status: str
    results: list[EvaluateReactionsResult]

class EvaluateReactionsResponse(BaseResponse):
    result: EvaluateReactionsResult


@register_wrapper(
    name="evaluate_reactions",
    input_class=EvaluateReactionsInput,
    output_class=EvaluateReactionsOutput,
    response_class=EvaluateReactionsResponse
)
class EvaluateReactionsWrapper(BaseWrapper):
    """Wrapper class for EvaluateReactionss"""
    prefixes = ["evaluate_reactions"]

    def call_sync(self, input: EvaluateReactionsInput) -> EvaluateReactionsResponse:
        output = super().call_raw(input=input)
        response = self.convert_output_to_response(output)
        return response

    async def call_async(self, input: EvaluateReactionsInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> EvaluateReactionsResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: EvaluateReactionsOutput
                                   ) -> EvaluateReactionsResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results[0]
        }
        if response["result"]:
            response = EvaluateReactionsResponse(**response)
            return response
        else:
            return None