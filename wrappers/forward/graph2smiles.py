from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel


class ForwardG2SInput(BaseModel):
    model_name: str
    smiles: list[str]


class ForwardG2SResult(BaseModel):
    products: list[str]
    scores: list[float]


class ForwardG2SOutput(BaseModel):
    __root__: list[ForwardG2SResult]


class ForwardG2SResponse(BaseResponse):
    result: list[ForwardG2SResult]


@register_wrapper(
    name="forward_graph2smiles",
    input_class=ForwardG2SInput,
    output_class=ForwardG2SOutput,
    response_class=ForwardG2SResponse
)
class ForwardG2SWrapper(BaseWrapper):
    """Wrapper class for Forward Prediction with Graph2SMILES"""
    prefixes = ["forward/graph2smiles"]

    def call_raw(self, input: ForwardG2SInput) -> ForwardG2SOutput:
        input_as_dict = input.dict()
        model_name = input_as_dict["model_name"]

        response = self.session_sync.post(
            f"{self.prediction_url}/{model_name}",
            json=input_as_dict,
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = ForwardG2SOutput(__root__=output)

        return output

    def call_sync(self, input: ForwardG2SInput) -> ForwardG2SResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: ForwardG2SInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ForwardG2SResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ForwardG2SOutput
                                   ) -> ForwardG2SResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.__root__
        }
        response = ForwardG2SResponse(**response)

        return response
