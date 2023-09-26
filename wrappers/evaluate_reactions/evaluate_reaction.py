from pydantic import BaseModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class EvaluateReactionInput(BaseModel):
    forward_backend: str = "wldn5"
    forward_model_name: str = "pistachio"
    context_n_conditions: int = 10

    rxn_smiles: str | None


class ForwardOutcome(BaseModel):
    outcome: str
    score: float


class ForwardTarget(BaseModel):
    smiles: str
    rank: int
    score: float


class EvaluateReactionResult(BaseModel):
    context: str
    top_product: ForwardOutcome
    number_of_outcomes: int
    target: ForwardTarget


class EvaluateReactionOutput(BaseModel):
    error: str
    status: str
    results: list[EvaluateReactionResult]


class EvaluateReactionResponse(BaseResponse):
    result: list[EvaluateReactionResult]


@register_wrapper(
    name="evaluate_reaction",
    input_class=EvaluateReactionInput,
    output_class=EvaluateReactionOutput,
    response_class=EvaluateReactionResponse
)
class EvaluateContextWrapper(BaseWrapper):
    """Wrapper class for evaluating context"""
    prefixes = ["evaluate_reaction"]

    def call_raw(self, input: EvaluateReactionInput) -> EvaluateReactionOutput:
        response = self.session_sync.post(
            f"{self.prediction_url}/evaluate_reaction",
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = EvaluateReactionOutput(**output)

        return output

    def call_sync(self, input: EvaluateReactionInput) -> EvaluateReactionResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: EvaluateReactionInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> EvaluateReactionResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: EvaluateReactionOutput
                                   ) -> EvaluateReactionResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = EvaluateReactionResponse(**response)

        return response
