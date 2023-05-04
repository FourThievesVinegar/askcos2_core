from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class AtomMapIndigoInput(BaseModel):
    smiles: list[str]


class AtomMapIndigoOutput(BaseModel):
    mapped_smiles: list[str]


@register_wrapper(
    name="atom_map_indigo",
    input_class=AtomMapIndigoInput,
    output_class=AtomMapIndigoOutput
)
class AtomMapIndigoWrapper(BaseWrapper):
    """Wrapper class for Indigo Atom Mapper"""
    DEFAULT_PREDICTION_URL = "http://atom_map_indigo:9919/"
    prefixes = ["atom_map/indigo"]
    methods_to_bind = {
        "get_config": ["GET"],
        "get_doc": ["GET"],
        "call_sync": ["POST"],
        "call_async": ["POST"],
        "retrieve": ["GET"]
    }

    def __init__(self, config: dict):
        super().__init__(config)

    def call_sync(self, input: AtomMapIndigoInput) -> AtomMapIndigoOutput:
        return super().call_sync(input=input)

    async def call_async(self, input: AtomMapIndigoInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapIndigoOutput | None:
        return await super().retrieve(task_id=task_id)
