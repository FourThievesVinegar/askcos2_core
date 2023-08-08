from wrappers import register_wrapper
from wrappers.base import BaseWrapper, BaseResponse
from pydantic import BaseModel
import json
from typing import List, Dict

class TreeOptimizerInput(BaseModel):
    targets: list
    graph: dict


class TreeOptimizerResult(BaseModel):
    tree: dict


class TreeOptimizerOutput(BaseModel):
    error: str
    status: str
    results: list[list[TreeOptimizerResult]]


class TreeOptimizerResponse(BaseResponse):
    result: List[TreeOptimizerResult]

@register_wrapper(
    name="tree_optimizer",
    input_class=TreeOptimizerInput,
    output_class=TreeOptimizerOutput,
    response_class=TreeOptimizerResponse
)
class TreeOptimizerWrapper(BaseWrapper):
    """Wrapper class for Tree Optimizer"""
    prefixes = ["tree_optimizer"]

    def call_raw(self, input: TreeOptimizerInput) -> TreeOptimizerOutput:
        response = self.session_sync.post(
            self.prediction_url,
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()

        return output

    def call_sync(self, input: TreeOptimizerInput
                  ) -> TreeOptimizerResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: TreeOptimizerInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> TreeOptimizerResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: TreeOptimizerOutput
                                   ) -> TreeOptimizerResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output["results"][0]
        }
        response = TreeOptimizerResponse(**response)

        return response




