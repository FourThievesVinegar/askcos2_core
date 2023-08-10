from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel


class DescriptorInput(BaseModel):
    smiles: str

class DescriptorResult(BaseModel):
    smiles: str
    partial_charge: list[float]
    fukui_neu: list[float]
    fukui_elec: list[float]
    NMR: list[float]
    bond_order: list[float]
    bond_length: list[float]


class DescriptorOutput(BaseModel):
    error: str
    status: str
    results: list[DescriptorResult]

class DescriptorResponse(BaseResponse):
    result: DescriptorResult


@register_wrapper(
    name="descriptors",
    input_class=DescriptorInput,
    output_class=DescriptorOutput,
    response_class=DescriptorResponse
)
class DescriptorWrapper(BaseWrapper):
    """Wrapper class for Descriptors"""
    prefixes = ["descriptors"]
    def call_sync(self, input: DescriptorInput):
        output = super().call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: DescriptorInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> DescriptorResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: DescriptorOutput
                                   ) -> DescriptorResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results[0]
        }
        if response["result"]:
            response = DescriptorResponse(**response)
            return response
        else:
            return None
