from pydantic import BaseModel, Field
from schemas.base import LowerCamelAliasModel
from typing import Literal
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class RetroATInput(LowerCamelAliasModel):
    model_name: Literal["pistachio_23Q2", "cas"] = "pistachio_23Q2"
    smiles: list[str]


class RetroATResult(BaseModel):
    reactants: list[str] = Field(alias="products")
    scores: list[float]


class RetroATOutput(BaseModel):
    __root__: list[RetroATResult]


class RetroATResponse(BaseResponse):
    result: list[RetroATResult]


@register_wrapper(
    name="retro_augmented_transformer",
    input_class=RetroATInput,
    output_class=RetroATOutput,
    response_class=RetroATResponse
)
class RetroATWrapper(BaseWrapper):
    """Wrapper class for Retro Prediction with Augmented Transformer"""
    prefixes = ["retro/augmented_transformer"]

    def call_raw(self, input: RetroATInput) -> RetroATOutput:
        input_as_dict = input.dict()
        model_name = input_as_dict["model_name"]

        response = self.session_sync.post(
            f"{self.prediction_url}/{model_name}",
            json=input_as_dict,
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = RetroATOutput(__root__=output)

        return output

    def call_sync(self, input: RetroATInput) -> RetroATResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: RetroATInput, priority: int = 0) -> str:
        from askcos2_celery.tasks import retro_task
        async_result = retro_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> RetroATResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: RetroATOutput
                                   ) -> RetroATResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.__root__
        }
        response = RetroATResponse(**response)

        return response
