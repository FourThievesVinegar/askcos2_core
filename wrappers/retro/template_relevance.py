from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel
from typing import Any


class RetroTemplRelInput(BaseModel):
    model_name: str
    smiles: list[str]


class RetroTemplRelResult(BaseModel):
    templates: list[dict[str, Any]]
    reactants: list[str]
    scores: list[float]


class RetroTemplRelOutput(BaseModel):
    __root__: list[RetroTemplRelResult]


@register_wrapper(
    name="retro_template_relevance",
    input_class=RetroTemplRelInput,
    output_class=RetroTemplRelOutput
)
class RetroTemplRelWrapper(BaseWrapper):
    """Wrapper class for One-step Retrosynthesis with Template Relevance Model"""
    prefixes = ["retro/template_relevance"]

    def call_sync(self, input: RetroTemplRelInput) -> RetroTemplRelOutput:
        input_as_dict = input.dict()
        model_name = input_as_dict["model_name"]

        response = self.session_sync.post(
            f"{self.prediction_url}/{model_name}",
            json=input_as_dict,
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = RetroTemplRelOutput(__root__=output)

        return output

    async def call_async(self, input: RetroTemplRelInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> RetroTemplRelOutput | None:
        return await super().retrieve(task_id=task_id)
