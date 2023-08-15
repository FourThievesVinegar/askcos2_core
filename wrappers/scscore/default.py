from pydantic import BaseModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class SCScoreInput(BaseModel):
    smiles: str


class SCScoreOutput(BaseModel):
    error: str
    status: str
    results: float


class SCScoreResponse(BaseResponse):
    result: float


@register_wrapper(
    name="scscore",
    input_class=SCScoreInput,
    output_class=SCScoreOutput,
    response_class=SCScoreResponse
)
class SCScoreWrapper(BaseWrapper):
    """Wrapper class for SCScore"""
    prefixes = ["scscore"]

    def call_sync(self, input: SCScoreInput) -> SCScoreResponse:
        output = super().call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: SCScoreInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> SCScoreResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: SCScoreOutput
                                   ) -> SCScoreResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = SCScoreResponse(**response)

        return response
