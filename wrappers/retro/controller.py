from pydantic import BaseModel, Field
from schemas.base import LowerCamelAliasModel
from scipy.special import softmax
from typing import Any, Literal
from utils.registry import get_util_registry
from wrappers import register_wrapper
from wrappers.base import BaseResponse, BaseWrapper
from wrappers.retro.augmented_transformer import RetroATInput, RetroATResponse
from wrappers.retro.graph2smiles import RetroG2SInput, RetroG2SResponse
from wrappers.retro.template_relevance import RetroTemplRelInput, RetroTemplRelResponse
from wrappers.registry import get_wrapper_registry


class AttributeFilter(BaseModel):
    name: str
    logic: Literal[">", ">=", "<", "<=", "=="]
    value: int | float


class RetroInput(LowerCamelAliasModel):
    backend: Literal[
        "augmented_transformer", "graph2smiles", "template_relevance"
    ] = Field(
        default="template_relevance",
        description="backend for one-step retrosynthesis"
    )
    model_name: str = Field(
        default="reaxys",
        description="backend model name for one-step retrosynthesis"
    )
    smiles: list[str] = Field(
        description="list of target SMILES",
        example=["CS(=N)(=O)Cc1cccc(Br)c1", "CN(C)CCOC(c1ccccc1)c1ccccc1"]
    )
    max_num_templates: int = Field(
        default=1000,
        description="number of templates to consider; "
                    "used for template_relevance only"
    )
    max_cum_prob: float = Field(
        default=0.995,
        description="maximum cumulative probability of templates; "
                    "used for template_relevance only"
    )
    attribute_filter: list[AttributeFilter] = Field(
        default_factory=list,
        description="template attribute filter to apply before template application; "
                    "used for template_relevance only",
        example=[]
    )


class RetroOutput(BaseModel):
    placeholder: str


class RetroResult(BaseModel):
    outcome: str
    model_score: float
    normalized_model_score: float
    template: dict[str, Any] | None


class RetroResponse(BaseResponse):
    result: list[list[RetroResult]]


@register_wrapper(
    name="retro_controller",
    input_class=RetroInput,
    output_class=RetroOutput,
    response_class=RetroResponse
)
class RetroController(BaseWrapper):
    """Retro Controller"""
    prefixes = ["retro/controller", "retro"]
    backend_wrapper_names = {
        "augmented_transformer": "retro_augmented_transformer",
        "graph2smiles": "retro_graph2smiles",
        "template_relevance": "retro_template_relevance"
    }

    def __init__(self):
        pass        # TODO: proper inheritance

    def call_sync(self, input: RetroInput) -> RetroResponse:
        """
        Endpoint for synchronous call to the retro controller,
        which dispatches the call to respective one-step retro backend service
        """
        cache_controller = get_util_registry().get_util(module="cache_controller")
        try:
            response = cache_controller.get(module_name=self.name, input=input)
            if isinstance(response, dict):
                response = RetroResponse(**response)
        except KeyError:
            module = self.backend_wrapper_names[input.backend]
            wrapper = get_wrapper_registry().get_wrapper(module=module)

            wrapper_input = self.convert_input(
                input=input, backend=input.backend)
            wrapper_response = wrapper.call_sync(wrapper_input)
            response = self.convert_response(
                wrapper_response=wrapper_response, backend=input.backend)

            cache_controller.add(module_name=self.name, input=input, response=response)

        return response

    async def call_async(self, input: RetroInput, priority: int = 0) -> str:
        """
        Endpoint for asynchronous call to the retro controller,
        which dispatches the call to respective one-step retro backend service
        """
        from askcos2_celery.tasks import retro_task
        async_result = retro_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> RetroResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_input(
        input: RetroInput, backend: str
    ) -> RetroATInput | RetroG2SInput | RetroTemplRelInput:
        if backend == "augmented_transformer":
            wrapper_input = RetroATInput(
                model_name=input.model_name,
                smiles=input.smiles
            )
        elif backend == "graph2smiles":
            wrapper_input = RetroG2SInput(
                model_name=input.model_name,
                smiles=input.smiles
            )
        elif backend == "template_relevance":
            wrapper_input = RetroTemplRelInput(
                model_name=input.model_name,
                smiles=input.smiles,
                max_num_templates=input.max_num_templates,
                max_cum_prob=input.max_cum_prob,
                attribute_filter=input.attribute_filter
            )
        else:
            raise ValueError(f"Unsupported retro backend: {backend}!")

        return wrapper_input

    @staticmethod
    def convert_response(
        wrapper_response:
            RetroATResponse | RetroG2SResponse | RetroTemplRelResponse,
        backend: str
    ) -> RetroResponse:
        status_code = wrapper_response.status_code
        message = wrapper_response.message
        if backend in ["augmented_transformer", "graph2smiles"]:
            # list[dict] -> list[list[dict]]
            result = []
            for result_per_smi in wrapper_response.result:
                normalized_scores = softmax(result_per_smi.scores)
                result.append(
                    [{
                        "outcome": outcome,
                        "model_score": score,
                        "normalized_model_score": float(normalized_score)
                    } for outcome, score, normalized_score in zip(
                        result_per_smi.reactants,
                        result_per_smi.scores,
                        normalized_scores
                    )]
                )
        elif backend == "template_relevance":
            # list[dict] -> list[list[dict]]
            result = []
            for result_per_smi in wrapper_response.result:
                # denominator = sum(result_per_smi.scores)  # we can re-normalize here?
                result.append(
                    [{
                        "outcome": outcome,
                        "model_score": score,
                        "normalized_model_score": score,
                        "template": template
                    } for outcome, score, template in zip(
                        result_per_smi.reactants,
                        result_per_smi.scores,
                        result_per_smi.templates
                    )]
                )
        else:
            raise ValueError(f"Unsupported retro backend: {backend}!")

        response = {
            "status_code": status_code,
            "message": message,
            "result": result
        }
        response = RetroResponse(**response)

        return response
