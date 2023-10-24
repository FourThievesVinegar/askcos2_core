from adapters import register_adapter
from pydantic import BaseModel, Field
from wrappers.registry import get_wrapper_registry
from wrappers.evaluate_reactions.evaluate_reaction import (
    EvaluateReactionInput,
    EvaluateReactionResponse
)


class V1EvaluationInput(BaseModel):
    smiles: str = Field(
        description="SMILES string of reaction to evaluate",
        example="CS(=N)(=O)Cc1cccc(Br)c1.Nc1cc(Br)ccn1>>CS(=N)(=O)Cc1cccc(Nc2cc(Br)ccn2)c1"
    )
    context_num_results: int = Field(
        default=10,
        description="number of context recommendations to evaluate",
        example=5
    )
    context_model_version: int = Field(
        default=1,
        description="which context model to use",
    )
    forward_num_results: int = Field(
        default=100,
        description="number of outcome predictions to search",
        example=20
    )
    forward_model_training_set: str = Field(
        default="uspto_500k",
        description="forward model training set to use"
    )
    forward_model_version: str = Field(
        default="",
        description="forward model version to use",
    )


class V1EvaluationAsyncReturn(BaseModel):
    request: V1EvaluationInput
    task_id: str


class V1Evaluation(BaseModel):
    context: list[str]
    top_product: dict
    number_of_outcomes: int
    target: dict


class V1EvaluationResult(BaseModel):
    rank: int
    smiles: str
    context: list[str]
    raw_context: list
    evaluation: V1Evaluation


class V1EvaluationResults(BaseModel):
    __root__: list[V1EvaluationResult]


@register_adapter(
    name="v1_evaluation",
    input_class=V1EvaluationInput,
    result_class=V1EvaluationResults
)
class V1EvaluationAdapter:
    """V1 Evaluation Adapter"""
    prefix = "legacy/evaluation"
    methods = ["POST"]

    def call_sync(self, input: V1EvaluationInput) -> V1EvaluationResults:
        wrapper = get_wrapper_registry().get_wrapper(
            module="evaluate_reaction"
        )
        wrapper_input = self.convert_input(input=input)
        wrapper_response = wrapper.call_sync(wrapper_input)
        response = self.convert_response(wrapper_response=wrapper_response)

        return response

    async def __call__(self, input: V1EvaluationInput, priority: int = 0
                       ) -> V1EvaluationAsyncReturn:
        from askcos2_celery.tasks import legacy_task

        async_result = legacy_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        async_return = V1EvaluationAsyncReturn(
            request=input, task_id=task_id)

        return async_return

    @staticmethod
    def convert_input(input: V1EvaluationInput) -> EvaluateReactionInput:
        wrapper_input = EvaluateReactionInput(
            rxn_smiles=input.smiles,
            forward_backend="wldn5",
            forward_model_name=input.forward_model_training_set,
            context_n_conditions=input.context_num_results
        )

        return wrapper_input

    @staticmethod
    def convert_response(wrapper_response: EvaluateReactionResponse
                         ) -> V1EvaluationResults:
        res = []
        for i, evaluation in enumerate(wrapper_response.result):
            evaluation.context = [evaluation.context]
            dict_evaluation = evaluation.dict()
            dict_evaluation["top_product"]["smiles"] = \
                dict_evaluation["top_product"].pop("outcome")

            v1_evaluation = V1EvaluationResult(
                rank=i+1,
                smiles="",
                context=evaluation.context,
                raw_context=[],
                evaluation=dict_evaluation
            )
            res.append(v1_evaluation)

        result = V1EvaluationResults(__root__=res)

        return result
