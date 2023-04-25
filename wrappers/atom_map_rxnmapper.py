from wrappers.base import BaseWrapper
from pydantic import BaseModel


class AtomMapRXNMapperInput(BaseModel):
    smiles: list[str]


class AtomMapRXNMapperOutput(BaseModel):
    mapped_smiles: list[str]


class AtomMapRXNMapperWrapper(BaseWrapper):
    DEFAULT_PREDICTION_URL = "http://0.0.0.0:7902/"
    DEFAULT_MANAGEMENT_URL = None
    input_class = AtomMapRXNMapperInput
    output_class = AtomMapRXNMapperOutput
    hparams_class = None

    name = "atom_map_rxnmapper"
    prefixes = ["/atom_map/rxnmapper"]
    methods_to_bind = ["call_sync", "call_async", "retrieve"]

    def __init__(self, config: dict, prediction_url: str = None,
                 management_url: str = None):
        super().__init__(config, prediction_url, management_url)

    def call_sync(self, input: AtomMapRXNMapperInput) -> AtomMapRXNMapperOutput:
        return super().call_sync(input=input)

    def call_async(self, input: AtomMapRXNMapperInput, priority: int
                   ) -> tuple[str, AtomMapRXNMapperInput]:
        return super().call_async(input=input, priority=priority)

    def retrieve(self, task_id: str) -> AtomMapRXNMapperOutput | None:
        return super().retrieve(task_id=task_id)
