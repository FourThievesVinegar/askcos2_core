from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel
from typing import Any, Literal


class AttributeFilter(BaseModel):
    name: str
    logic: Literal[">", ">=", "<", "<=", "=="]
    value: int | float


class RetroTemplRelInput(BaseModel):
    model_name: str
    smiles: list[str]
    max_num_templates: int = 1000
    max_cum_prob: float = 0.999
    attribute_filter: list[AttributeFilter] = []


class RetroTemplRelResult(BaseModel):
    templates: list[dict[str, Any]]
    reactants: list[str]
    scores: list[float]


class RetroTemplRelOutput(BaseModel):
    __root__: list[RetroTemplRelResult]


class RetroTemplRelResponse(BaseResponse):
    result: list[RetroTemplRelResult]


@register_wrapper(
    name="retro_template_relevance",
    input_class=RetroTemplRelInput,
    output_class=RetroTemplRelOutput,
    response_class=RetroTemplRelResponse
)
class RetroTemplRelWrapper(BaseWrapper):
    """Wrapper class for One-step Retrosynthesis with Template Relevance Model"""
    prefixes = ["retro/template_relevance"]

    def call_raw(self, input: RetroTemplRelInput) -> RetroTemplRelOutput:
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

    def call_sync(self, input: RetroTemplRelInput) -> RetroTemplRelResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: RetroTemplRelInput, priority: int = 0) -> str:
        from askcos2_celery.tasks import retro_task
        async_result = retro_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> RetroTemplRelResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: RetroTemplRelOutput
                                   ) -> RetroTemplRelResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.__root__
        }
        response = RetroTemplRelResponse(**response)

        return response
