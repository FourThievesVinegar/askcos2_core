from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class ForwardATInput(BaseModel):
    model_name: str
    smiles: list[str]


class ForwardATResult(BaseModel):
    products: list[str]
    scores: list[float]


class ForwardATOutput(BaseModel):
    __root__: list[ForwardATResult]


@register_wrapper(
    name="forward_augmented_transformer",
    input_class=ForwardATInput,
    output_class=ForwardATOutput
)
class ForwardATWrapper(BaseWrapper):
    """Wrapper class for Forward Prediction with Augmented Transformer"""
    prefixes = ["forward/augmented_transformer"]

    def call_sync(self, input: ForwardATInput) -> ForwardATOutput:
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

    async def call_async(self, input: ForwardATInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ForwardATOutput | None:
        return await super().retrieve(task_id=task_id)
