from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class SiteSelectivityInput(BaseModel):
    smiles: list[str]


class SiteSelectivityResult(BaseModel):
    smiles: str
    index: int
    task: str
    atom_scores: list[float]


class SiteSelectivityOutput(BaseModel):
    error: str
    status: str
    results: list[SiteSelectivityResult]


@register_wrapper(
    name="site_selectivity",
    input_class=SiteSelectivityInput,
    output_class=SiteSelectivityOutput
)
class SiteSelectivityWrapper(BaseWrapper):
    """Wrapper class for Site Selectivity"""
    prefixes = ["site_selectivity"]

    def call_sync(self, input: SiteSelectivityInput) -> SiteSelectivityOutput:
        return super().call_sync(input=input)

    async def call_async(self, input: SiteSelectivityInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> SiteSelectivityOutput | None:
        return await super().retrieve(task_id=task_id)
