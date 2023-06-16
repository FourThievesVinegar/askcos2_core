from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class FastFilterInput(BaseModel):
    smiles: list[str]


class FastFilterOutcome(BaseModel):
    smiles: str
    template_ids: list
    num_examples: int


class FastFilterResult(BaseModel):
    rank: float
    outcome: FastFilterOutcome
    score: float
    prob: float


class FastFilterOutput(BaseModel):
    error: str
    status: str
    results: list[list[list[FastFilterResult]]]


@register_wrapper(
    name="fast_filter",
    input_class=FastFilterInput,
    output_class=FastFilterOutput
)
class FastFilterWrapper(BaseWrapper):
    """Wrapper class for Fast Filter"""
    prefixes = ["fast_filter"]

    def call_sync(self, input: FastFilterInput) -> FastFilterOutput:
        self.prediction_url = f"{self.prediction_url}/fast_filter_evaluate"

        return super().call_sync(input=input)

    async def call_async(self, input: FastFilterInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> FastFilterOutput | None:
        return await super().retrieve(task_id=task_id)
