from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel
from typing import Literal


class GeneralSelectivityInput(BaseModel):
    smiles: str
    atom_map_backend: Literal["indigo", "rxnmapper", "wln"] = "wln"


class GeneralSelectivityResult(BaseModel):
    smiles: str
    prob: float
    rank: int


class GeneralSelectivityOutput(BaseModel):
    error: str
    status: str
    results: list[GeneralSelectivityResult]


class GeneralSelectivityResponse(BaseResponse):
    result: list[GeneralSelectivityResult]


@register_wrapper(
    name="general_selectivity_gnn",
    input_class=GeneralSelectivityInput,
    output_class=GeneralSelectivityOutput,
    response_class=GeneralSelectivityResponse
)
class GeneralSelectivityGNNWrapper(BaseWrapper):
    """Wrapper class for general selectivity, GNN model"""
    prefixes = ["general_selectivity/gnn"]

    def call_raw(self, input: GeneralSelectivityInput) -> GeneralSelectivityOutput:
        response = self.session_sync.post(
            f"{self.prediction_url}/general_selectivity",
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = GeneralSelectivityOutput(**output)

        return output

    def call_sync(self, input: GeneralSelectivityInput) -> GeneralSelectivityResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: GeneralSelectivityInput, priority: int = 0
                         ) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> GeneralSelectivityResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: GeneralSelectivityOutput
                                   ) -> GeneralSelectivityResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = GeneralSelectivityResponse(**response)

        return response
