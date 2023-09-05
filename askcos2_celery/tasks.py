from celery import shared_task
from adapters.registry import get_adapter_registry
from wrappers.registry import get_wrapper_registry


@shared_task
def base_task(module: str, input: dict) -> dict:
    """Celery tasks must have json serializable inputs/outputs"""
    wrapper = get_wrapper_registry().get_wrapper(module=module)

    # Reconstruct Input object from, and convert Output object to dict
    input = wrapper.input_class(**input)
    response = wrapper.call_sync(input).dict()

    return response


@shared_task
def legacy_task(module: str, input: dict) -> dict:
    """Celery tasks must have json serializable inputs/outputs"""
    adapter = get_adapter_registry().get_adapter(module=module)

    # Reconstruct Input object from, and convert Output object to dict
    input = adapter.input_class(**input)
    response = adapter.call_sync(input).dict()

    return response
