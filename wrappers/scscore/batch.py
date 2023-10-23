from pydantic import Field
from schemas.base import LowerCamelAliasModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from wrappers.scscore.default import SCScoreInput, SCScoreOutput


class SCScoreBatchInput(LowerCamelAliasModel):
    smiles: list[str] = Field(
        description="list of SMILES for SCScore calculation. "
                    "Can be multi-molecule (separated by .), in which case "
                    "the max SCScore among the molecules will be returned",
        example=["c1ccccc1", "CCCBr"]
    )


class SCScoreBatchResponse(BaseResponse):
    result: dict[str, float] | None


@register_wrapper(
    name="scscore_batch",
    input_class=SCScoreBatchInput,
    output_class=SCScoreOutput,
    response_class=SCScoreBatchResponse
)
class SCScoreWrapper(BaseWrapper):
    """Wrapper class for SCScore"""
    prefixes = ["scscore/batch"]

    def call_sync(self, input: SCScoreBatchInput) -> SCScoreBatchResponse:
        """
        Endpoint for synchronous call to batch SCScorer.
        https://pubs.acs.org/doi/10.1021/acs.jcim.7b00622
        """
        results = {}
        for smi in input.smiles:
            output = self.call_raw(input=SCScoreInput(smiles=smi))
            score = output.results
            results[smi] = score

        response = self.convert_results_to_response(results)

        return response

    async def call_async(self, input: SCScoreBatchInput, priority: int = 0) -> str:
        """
        Endpoint for asynchronous call to batch SCScorer.
        https://pubs.acs.org/doi/10.1021/acs.jcim.7b00622
        """
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> SCScoreBatchResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_results_to_response(results: dict[str, float]
                                    ) -> SCScoreBatchResponse:

        response = SCScoreBatchResponse(
            status_code=200,
            message="",
            result=results
        )

        return response
