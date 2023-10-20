from pydantic import BaseModel, Field
from schemas.base import LowerCamelAliasModel
from wrappers import register_wrapper
from wrappers.base import BaseWrapper


class LegacyDescriptorsInput(LowerCamelAliasModel):
    smiles: str = Field(
        description="SMILES of the molecule",
        example="CCC"
    )


class LegacyDescriptorsOutput(BaseModel):
    smiles: str
    partial_charge: list[float]
    fukui_neu: list[float]
    fukui_elec: list[float]
    NMR: list[float]
    bond_order: list[float]
    bond_length: list[float]


# trivial; there doesn't seem to be postprocessing for descriptors
LegacyDescriptorsResponse = LegacyDescriptorsOutput


@register_wrapper(
    name="legacy_descriptors",
    input_class=LegacyDescriptorsInput,
    output_class=LegacyDescriptorsOutput,
    response_class=LegacyDescriptorsResponse
)
class LegacyDescriptorsWrapper(BaseWrapper):
    """Wrapper class for Legacy Descriptors"""
    prefixes = ["legacy/descriptors"]

    def call_raw(self, input: LegacyDescriptorsInput) -> LegacyDescriptorsOutput:
        response = self.session_sync.post(
            self.prediction_url,
            json=[input.smiles],
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = LegacyDescriptorsOutput(**output)

        return output

    def call_sync(self, input: LegacyDescriptorsInput) -> LegacyDescriptorsResponse:
        """
        Endpoint for synchronous call to QM descriptors service
        """
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: LegacyDescriptorsInput, priority: int = 0
                         ) -> str:
        """
        Endpoint for asynchronous call to QM descriptors service
        """
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> LegacyDescriptorsResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: LegacyDescriptorsOutput
                                   ) -> LegacyDescriptorsResponse:
        # trivial; there doesn't seem to be postprocessing for descriptors
        return output
