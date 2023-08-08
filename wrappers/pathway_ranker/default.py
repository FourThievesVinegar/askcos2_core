from wrappers import register_wrapper
from wrappers.base import BaseWrapper, BaseResponse
from fastapi.encoders import jsonable_encoder
from pydantic import BaseModel
import json
import requests
from typing import List, Dict

class PathwayRankerInput(BaseModel):
    tree: List[Dict]


class PathwayRankerResult(BaseModel):
    scores: list
    encoded_trees: list
    clusters: list

class PathwayRankerOutput(BaseModel):
    error: str
    status: str
    results: list[PathwayRankerResult]

class PathwayRankerResponse(BaseResponse):
    result: List[PathwayRankerResult]

@register_wrapper(
    name="pathway_ranker",
    input_class=PathwayRankerInput,
    output_class=PathwayRankerOutput,
    response_class=PathwayRankerResponse
)
class PathwayRankerWrapper(BaseWrapper):
    """Wrapper class for Pathway Ranker"""
    prefixes = ["pathway_ranker"]
    def call_sync(self, input: PathwayRankerInput) -> PathwayRankerResponse:
        output = super().call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: PathwayRankerInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> PathwayRankerResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: PathwayRankerOutput
                                   ) -> PathwayRankerResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = PathwayRankerResponse(**response)

        return response