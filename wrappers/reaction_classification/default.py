from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class ReactionClassificationInput(BaseModel):
    smiles: list[str]


class ReactionClassificationResult(BaseModel):
    rank: float
    reaction_num: str
    reaction_name: str
    reaction_classnum: float
    reaction_classname: str
    reaction_superclassnum: float
    reaction_superclassname: str
    prediction_certainty: float


class ReactionClassificationOutput(BaseModel):
    error: str
    status: str
    results: list[ReactionClassificationResult]


@register_wrapper(
    name="reaction_classification",
    input_class=ReactionClassificationInput,
    output_class=ReactionClassificationOutput
)
class ReactionClassificationWrapper(BaseWrapper):
    """Wrapper class for Reaction Classification"""
    prefixes = ["reaction_classification"]

    def call_sync(self, input: ReactionClassificationInput) -> ReactionClassificationOutput:
        return super().call_sync(input=input)

    async def call_async(self, input: ReactionClassificationInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ReactionClassificationOutput | None:
        return await super().retrieve(task_id=task_id)
