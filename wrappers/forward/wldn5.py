import copy
from pydantic import BaseModel
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper


class ForwardWLDN5Input(BaseModel):
    model_name: str
    reactants: str | list[str]
    contexts: list[str] | None = None


class ForwardWLDN5Outcome(BaseModel):
    smiles: str
    template_ids: list = None
    num_examples: int = 0


class ForwardWLDN5Result(BaseModel):
    rank: int
    outcome: ForwardWLDN5Outcome
    score: float
    prob: float
    mol_wt: float


class ForwardWLDN5Output(BaseModel):
    status: str
    error: str
    results: list[list[ForwardWLDN5Result]]


class ForwardWLDN5Response(BaseResponse):
    result: list[list[ForwardWLDN5Result]]


@register_wrapper(
    name="forward_wldn5",
    input_class=ForwardWLDN5Input,
    output_class=ForwardWLDN5Output,
    response_class=ForwardWLDN5Response
)
class ForwardWLDN5Wrapper(BaseWrapper):
    """Wrapper class for Forward Prediction with WLDN5"""
    prefixes = ["forward/wldn5"]

    def call_raw(self, input: ForwardWLDN5Input) -> ForwardWLDN5Output:
        if isinstance(input.reactants, str):
            input.reactants = [input.reactants]

        input_as_dict = copy.deepcopy(input.dict())
        results = []
        for i, reactants in enumerate(input.reactants):
            input_as_dict["reactants"] = reactants
            response = self.session_sync.post(
                f"{self.prediction_url}",
                json=input_as_dict,
                timeout=self.config["deployment"]["timeout"]
            )
            result = response.json()["results"]
            results.extend(result)

        output = response.json()            # hardcoding using the last response status
        output["results"] = results
        output = ForwardWLDN5Output(**output)

        return output

    def call_sync(self, input: ForwardWLDN5Input) -> ForwardWLDN5Response:
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: ForwardWLDN5Input, priority: int = 0) -> str:
        from askcos2_celery.tasks import forward_task
        async_result = forward_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> ForwardWLDN5Response | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: ForwardWLDN5Output
                                   ) -> ForwardWLDN5Response:
        response = {
            "status_code": 200,
            "message": "",
            "result": output.results
        }
        response = ForwardWLDN5Response(**response)

        return response
