from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel
from typing import Optional


class ContextRecommenderInput(BaseModel):
    smiles: str
    n_conditions: Optional[int] = 10


# class ContextRecommenderResult(BaseModel):
#     agents: list
#     temperature: float
#     score: float


class ContextRecommenderOutput(BaseModel):
    __root__: list[dict]
    # __root__: list[ContextRecommenderResult]



# FP model condition
@register_wrapper(
    name="context_recommender_pr_fp", # So this must be "context_recommender", not context_recommender_fp? But why? (Aaron Chen)
    input_class=ContextRecommenderInput,
    output_class=ContextRecommenderOutput
)
class ContextRecommenderWrapper(BaseWrapper):
    """Wrapper class for Context Recommender Predict FP"""
    prefixes = ["context_recommender/v2/predict/FP"]

    def call_sync(self, input: ContextRecommenderInput):
        response = self.session_sync.post(
            f"{self.prediction_url}/api/v2/predict/FP",
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = ContextRecommenderOutput(__root__=output)
        return output

    async def call_async(self, input: ContextRecommenderInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ContextRecommenderOutput | None:
        return await super().retrieve(task_id=task_id)


