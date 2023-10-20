from pydantic import BaseModel, Field
from schemas.base import LowerCamelAliasModel
from typing import Literal
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class EvaluateReactionInput(LowerCamelAliasModel):
    forward_backend: Literal["wldn5", "augmented_transformer", "graph2smiles"] = Field(
        default="wldn5",
        description="Backend for forward prediction"
    )
    forward_model_name: str = Field(
        default="pistachio",
        description="Backend model name for forward prediction"
    )
    context_n_conditions: int = Field(
        default=10,
        description="Number of conditions to predict from the context recommender"
    )

    rxn_smiles: str = Field(
        description="Reaction SMILES to be evaluated (w/o reagents)",
        example="CS(=N)(=O)Cc1cccc(Br)c1.Nc1cc(Br)ccn1>>"
                "CS(=N)(=O)Cc1cccc(Nc2cc(Br)ccn2)c1"
    )


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
    result: list[EvaluateReactionResult] | None


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
        """
        Endpoint for synchronous call to reaction evaluation workflow.
        Predict reaction conditions and then run forward evaluation for each condition.
        """
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: EvaluateReactionInput, priority: int = 0) -> str:
        """
        Endpoint for asynchronous call to reaction evaluation workflow.
        Predict reaction conditions and then run forward evaluation for each condition.
        """
        from askcos2_celery.tasks import evaluate_reactions_task
        async_result = evaluate_reactions_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> EvaluateReactionResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: EvaluateReactionOutput
                                   ) -> EvaluateReactionResponse:
        if output.status == "SUCCESS":
            status_code = 200
            message = ""
            result = output.results
        else:
            status_code = 500
            message = f"Backend error encountered in evaluate_reaction " \
                      f"with the following error message {output.error}"
            result = None

        response = EvaluateReactionResponse(
            status_code=status_code,
            message=message,
            result=result
        )

        return response
