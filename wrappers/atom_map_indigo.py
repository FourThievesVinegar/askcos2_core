from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class AtomMapIndigoInput(BaseModel):
    smiles: list[str]


class AtomMapIndigoOutput(BaseModel):
    mapped_smiles: list[str]


@register_wrapper(name="atom_map_indigo")
class AtomMapIndigoWrapper(BaseWrapper):
    DEFAULT_PREDICTION_URL = "http://0.0.0.0:7901/"
    DEFAULT_MANAGEMENT_URL = None
    input_class = AtomMapIndigoInput
    output_class = AtomMapIndigoOutput

    name = "atom_map_indigo"
    prefixes = ["atom_map/indigo"]
    methods_to_bind = {
        "get_metadata": ["GET"],
        "call_sync": ["POST"],
        "call_async": ["POST"],
        "retrieve": ["GET"]
    }

    def __init__(self, config: dict, prediction_url: str = None,
                 management_url: str = None):
        super().__init__(config, prediction_url, management_url)

    def call_sync(self, input: input_class) -> output_class:
        return super().call_sync(input=input)

    async def call_async(self, input: input_class, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapIndigoOutput | None:
        return await super().retrieve(task_id=task_id)
