from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel
from typing import Optional


class ContextRecommenderInput(BaseModel):
    smiles: str
    n_conditions: Optional[int] = 10


class ContextRecommenderOutput(BaseModel):
    __root__: list[dict]


class ContextRecommenderResponse(BaseResponse):
    result: list[dict]


# FP model condition
@register_wrapper(
    name="context_recommender_pr_graph",
    input_class=ContextRecommenderInput,
    output_class=ContextRecommenderOutput,
    response_class=ContextRecommenderResponse
)
class ContextRecommenderWrapper(BaseWrapper):
    """Wrapper class for Context Recommender Predict GRAPH"""
    prefixes = ["context_recommender/v2/predict/GRAPH"]

    def call_raw(self, input: ContextRecommenderInput) -> ContextRecommenderOutput:
        response = self.session_sync.post(
            f"{self.prediction_url}/api/v2/predict/GRAPH",
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = ContextRecommenderOutput(__root__=output)

        return output

    def call_sync(self, input: ContextRecommenderInput) -> ContextRecommenderResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: ContextRecommenderInput, priority: int = 0
                         ) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ContextRecommenderOutput | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ContextRecommenderOutput
                                   ) -> ContextRecommenderResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result_format": "json",
            "result": output.__root__
        }
        response = ContextRecommenderResponse(**response)

        return response
