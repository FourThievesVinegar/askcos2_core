from pydantic import BaseModel
from schemas.base import LowerCamelAliasModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class EvaluateContextInput(LowerCamelAliasModel):
    forward_backend: str = "wldn5"
    forward_model_name: str = "pistachio"

    reactant_smiles: str | None
    product_smiles: str | None
    rxn_smiles: str | None
    contexts: list[str]
    top_n: int = 10
    return_all_outcomes: bool = False
    isomeric_smiles: bool = True


class ForwardOutcome(BaseModel):
    outcome: str
    score: float


class ForwardTarget(BaseModel):
    smiles: str
    rank: int
    score: float


class EvaluateContextResult(BaseModel):
    context: str
    top_product: ForwardOutcome
    number_of_outcomes: int
    target: ForwardTarget


class EvaluateContextOutput(BaseModel):
    error: str
    status: str
    results: list[EvaluateContextResult]


class EvaluateContextResponse(BaseResponse):
    result: list[EvaluateContextResult]


@register_wrapper(
    name="evaluate_context",
    input_class=EvaluateContextInput,
    output_class=EvaluateContextOutput,
    response_class=EvaluateContextResponse
)
class EvaluateContextWrapper(BaseWrapper):
    """Wrapper class for evaluating context"""
    prefixes = ["evaluate_context"]

    def call_raw(self, input: EvaluateContextInput) -> EvaluateContextOutput:
        response = self.session_sync.post(
            f"{self.prediction_url}/evaluate_context",
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = EvaluateContextOutput(**output)

        return output

    def call_sync(self, input: EvaluateContextInput) -> EvaluateContextResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: EvaluateContextInput, priority: int = 0) -> str:
        from askcos2_celery.tasks import evaluate_reactions_task
        async_result = evaluate_reactions_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> EvaluateContextResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: EvaluateContextOutput
                                   ) -> EvaluateContextResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = EvaluateContextResponse(**response)

        return response
