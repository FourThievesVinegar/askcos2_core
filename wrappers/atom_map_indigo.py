from wrappers.base import TorchServeWrapper
from pydantic import BaseModel


class AtomMapIndigoInput(BaseModel):
    smiles: list[str]


class AtomMapIndigoOutput(BaseModel):
    mapped_smiles: list[str]


class AtomMapIndigoWrapper(TorchServeWrapper):
    DEFAULT_PREDICTION_URL = "http://0.0.0.0:7901/"
    DEFAULT_MANAGEMENT_URL = None
    input_class = AtomMapIndigoInput
    output_class = AtomMapIndigoOutput
    hparams_class = None

    name = "atom_map_indigo"
    prefixes = ["/atom_map/indigo"]
    methods_to_bind = ["call_sync", "call_async", "retrieve"]

    def __init__(self, config: dict, prediction_url: str = None,
                 management_url: str = None):
        super().__init__(config, prediction_url, management_url)

    def call_sync(self, input: AtomMapIndigoInput) -> AtomMapIndigoOutput:
        return super().call_sync(input=input)

    def call_async(self, input: AtomMapIndigoInput, priority: int
                   ) -> tuple[str, AtomMapIndigoInput]:
        return super().call_async(input=input, priority=priority)

    def retrieve(self, task_id: str) -> AtomMapIndigoOutput | None:
        return super().retrieve(task_id=task_id)
