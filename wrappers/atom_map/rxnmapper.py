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
    DEFAULT_PREDICTION_URL = "http://atom_map_rxnmapper_service:9918/rxnmapper"
    prefixes = ["atom_map/rxnmapper"]

    def call_sync(self, input: AtomMapRXNMapperInput) -> AtomMapRXNMapperOutput:
        return super().call_sync(input=input)

    async def call_async(self, input: AtomMapRXNMapperInput, priority: int = 0) -> str:
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> AtomMapRXNMapperOutput | None:
        return await super().retrieve(task_id=task_id)
