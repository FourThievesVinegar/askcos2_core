from pydantic import BaseModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class ForwardATInput(BaseModel):
    model_name: str
    smiles: list[str]


class ForwardATResult(BaseModel):
    products: list[str]
    scores: list[float]


class ForwardATOutput(BaseModel):
    __root__: list[ForwardATResult]


class ForwardATResponse(BaseResponse):
    result: list[ForwardATResult]


@register_wrapper(
    name="forward_augmented_transformer",
    input_class=ForwardATInput,
    output_class=ForwardATOutput,
    response_class=ForwardATResponse
)
class ForwardATWrapper(BaseWrapper):
    """Wrapper class for Forward Prediction with Augmented Transformer"""
    prefixes = ["forward/augmented_transformer"]

    def call_raw(self, input: ForwardATInput) -> ForwardATOutput:
        input_as_dict = input.dict()
        model_name = input_as_dict["model_name"]

        response = self.session_sync.post(
            f"{self.prediction_url}/{model_name}",
            json=input_as_dict,
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = ForwardATOutput(__root__=output)

        return output

    def call_sync(self, input: ForwardATInput) -> ForwardATResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: ForwardATInput, priority: int = 0) -> str:
        from askcos2_celery.tasks import forward_task
        async_result = forward_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> ForwardATResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ForwardATOutput
                                   ) -> ForwardATResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.__root__
        }
        response = ForwardATResponse(**response)

        return response
