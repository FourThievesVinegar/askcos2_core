from adapters import register_adapter
from pydantic import BaseModel
from wrappers.registry import get_wrapper_registry
from wrappers.descriptors.default import (
    DescriptorsInput,
    DescriptorsResult,
    DescriptorsResponse
)


class V1DescriptorsInput(BaseModel):
    smiles: str


class V1DescriptorsAsyncReturn(BaseModel):
    request: V1DescriptorsInput
    task_id: str


class V1DescriptorsResult(BaseModel):
    __root__: list[DescriptorsResult]


@register_adapter(
    name="v1_descriptors",
    input_class=V1DescriptorsInput,
    result_class=V1DescriptorsResult
)
class V1DescriptorsAdapter:
    """V1 Descriptors Adapter"""
    prefix = "legacy/descriptors"
    methods = ["POST"]

    def call_sync(self, input: V1DescriptorsInput) -> V1DescriptorsResult:
        wrapper = get_wrapper_registry().get_wrapper(module="descriptors")
        wrapper_input = self.convert_input(input=input)
        wrapper_response = wrapper.call_sync(wrapper_input)
        response = self.convert_response(wrapper_response=wrapper_response)

        return response

    async def __call__(self, input: V1DescriptorsInput, priority: int = 0
                       ) -> V1DescriptorsAsyncReturn:
        from askcos2_celery.tasks import legacy_task

        async_result = legacy_task.apply_async(
            args=(self.name, input.dict()), priority=priority)
        task_id = async_result.id

        async_return = V1DescriptorsAsyncReturn(
            request=input, task_id=task_id)

        return async_return

    @staticmethod
    def convert_input(input: V1DescriptorsInput) -> DescriptorsInput:
        wrapper_input = DescriptorsInput(smiles=[input.smiles])

        return wrapper_input

    @staticmethod
    def convert_response(wrapper_response: DescriptorsResponse
                         ) -> V1DescriptorsResult:
        result = wrapper_response.result[0]
        result = V1DescriptorsResult(__root__=result)

        return result
