from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from pydantic import BaseModel


class RXNMapperResult(BaseModel):
    mapped_rxn: str
    confidence: float


class AtomMapRXNMapperInput(BaseModel):
    smiles: list[str]


class AtomMapRXNMapperOutput(BaseModel):
    error: str
    status: str
    results: list[RXNMapperResult]


@register_wrapper(
    name="atom_map_rxnmapper",
    input_class=AtomMapRXNMapperInput,
    output_class=AtomMapRXNMapperOutput
)
class AtomMapRXNMapperWrapper(BaseWrapper):
    """Wrapper class for IBM RXNMapper"""
    DEFAULT_PREDICTION_URL = "http://atom_map_rxnmapper:9918/"
    prefixes = ["atom_map/rxnmapper"]
    methods_to_bind = {
        "get_config": ["GET"],
        "get_doc": ["GET"],
        "call_sync": ["POST"],
        "call_async": ["POST"],
        "retrieve": ["GET"]
    }

    def __init__(self, config: dict):
        super().__init__(config)

    def call_sync(self, input: AtomMapRXNMapperInput) -> AtomMapRXNMapperOutput:
        return super().call_sync(input=input)

    async def call_async(self, input: AtomMapRXNMapperInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapRXNMapperOutput | None:
        return await super().retrieve(task_id=task_id)
