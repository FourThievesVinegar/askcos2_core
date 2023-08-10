from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel
import json


class GeneralSelectivityInput(BaseModel):
    smiles: list[str]


class GeneralSelectivityResult(BaseModel):
    smiles: str
    prob: float
    rank: int


class GeneralSelectivityOutput(BaseModel):
    __root__: list[list[GeneralSelectivityResult]]

class GeneralSelectivityResponse(BaseResponse):
    result: list[list[GeneralSelectivityResult]]

@register_wrapper(
    name="general_selectivity_qm_no_reagent",
    input_class=GeneralSelectivityInput,
    output_class=GeneralSelectivityOutput,
    response_class=GeneralSelectivityResponse
)
class GeneralSelectivityWrapper(BaseWrapper):
    """Wrapper class for General Selectivity Qm No Reagent"""
    prefixes = ["general_selectivity/qm_no_reagent"]
    def call_raw(self, input: GeneralSelectivityInput)->GeneralSelectivityOutput:
        response = self.session_sync.post(
            f"{self.prediction_url}/qm_no_reagent",
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = GeneralSelectivityOutput(__root__=output["results"])
        return output

    def call_sync(self, input: GeneralSelectivityInput) -> GeneralSelectivityResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)
        return response

    async def call_async(self, input: GeneralSelectivityInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> GeneralSelectivityResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: GeneralSelectivityOutput
                                   ) -> GeneralSelectivityResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.__root__
        }
        if response["result"]:
            response = GeneralSelectivityResponse(**response)
            return response
        else:
            return None



