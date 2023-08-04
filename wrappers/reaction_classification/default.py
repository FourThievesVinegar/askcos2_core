from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from pydantic import BaseModel


class ReactionClassificationInput(BaseModel):
    smiles: list[str]
    num_results: int = 10


class ReactionClassificationResult(BaseModel):
    rank: float
    reaction_num: str
    reaction_name: str
    reaction_classnum: float
    reaction_classname: str
    reaction_superclassnum: float
    reaction_superclassname: str
    prediction_certainty: float


class ReactionClassificationOutput(BaseModel):
    error: str
    status: str
    results: list[ReactionClassificationResult]


class ReactionClassificationResponse(BaseResponse):
    result: list[ReactionClassificationResult]


@register_wrapper(
    name="reaction_classification",
    input_class=ReactionClassificationInput,
    output_class=ReactionClassificationOutput,
    response_class=ReactionClassificationResponse
)
class ReactionClassificationWrapper(BaseWrapper):
    """Wrapper class for Reaction Classification"""
    prefixes = ["reaction_classification"]

    def call_raw(self, input: ReactionClassificationInput
                 ) -> ReactionClassificationOutput:
        response = self.session_sync.post(
            f"{self.prediction_url}/reaction_class",
            json=input.dict(),
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = ReactionClassificationOutput(**output)

        return output

    def call_sync(self, input: ReactionClassificationInput
                  ) -> ReactionClassificationResponse:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: ReactionClassificationInput, priority: int = 0
                         ) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> ReactionClassificationResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ReactionClassificationOutput
                                   ) -> ReactionClassificationResponse:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = ReactionClassificationResponse(**response)

        return response
