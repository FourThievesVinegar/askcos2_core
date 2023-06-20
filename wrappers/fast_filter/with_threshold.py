from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class FastFilterWThresholdInput(BaseModel):
    smiles: list[str]
    threshold: float


class FastFilterWThresholdResult(BaseModel):
    flag: bool
    score: float


class FastFilterWThresholdOutput(BaseModel):
    error: str
    status: str
    results: list[FastFilterWThresholdResult]


@register_wrapper(
    name="fast_filter_with_threshold",
    input_class=FastFilterWThresholdInput,
    output_class=FastFilterWThresholdOutput
)
class FastFilterWThresholdWrapper(BaseWrapper):
    """Wrapper class for Fast Filter with Threshold"""
    prefixes = ["fast_filter/with_threshold"]

    def call_sync(self, input: FastFilterWThresholdInput
                  ) -> FastFilterWThresholdOutput:
        self.prediction_url = f"{self.prediction_url}/filter_with_threshold"

        return super().call_sync(input=input)

    async def call_async(self, input: FastFilterWThresholdInput, priority: int = 0
                         ) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> FastFilterWThresholdOutput | None:
        return await super().retrieve(task_id=task_id)
