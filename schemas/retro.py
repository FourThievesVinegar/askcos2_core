from pydantic import Field, validator
from schemas.base import LowerCamelAliasModel
from typing import Any, Literal


class RetroBackendOption(LowerCamelAliasModel):
    retro_backend: Literal[
        "augmented_transformer", "graph2smiles", "template_relevance", "retrosim"
    ] = Field(
        default="template_relevance",
        description="backend for one-step retrosynthesis"
    )
    retro_model_name: str = Field(
        default="reaxys",
        description="backend model name for one-step retrosynthesis"
    )

    max_num_templates: int = Field(
        default=1000,
        description="number of templates to consider"
    )
    max_cum_prob: float = Field(
        default=0.995,
        description="maximum cumulative probability of templates"
    )
    attribute_filter: list[dict[str, Any]] = Field(
        default_factory=list,
        description="template attribute filter to apply before template application",
        example=[]
    )

    @validator("retro_model_name")
    def check_retro_model_name(cls, v, values):
        if "retro_backend" not in values:
            raise ValueError("retro_backend not supplied!")
        if values["retro_backend"] == "template_relevance":
            if v not in [
                "bkms_metabolic",
                "cas",
                "pistachio",
                "pistachio_ringbreaker",
                "reaxys",
                "reaxys_biocatalysis"
            ]:
                raise ValueError(
                    f"Unsupported retro_model_name {v} for template_relevance")
        elif values["retro_backend"] == "augmented_transformer":
            if v not in [
                "cas",
                "pistachio_23Q3",
                "USPTO_FULL"
            ]:
                raise ValueError(
                    f"Unsupported retro_model_name {v} for augmented_transformer")
        elif values["retro_backend"] == "graph2smiles":
            if v not in [
                "cas",
                "pistachio_23Q3",
                "USPTO_FULL"
            ]:
                raise ValueError(
                    f"Unsupported retro_model_name {v} for graph2smiles")
        elif values["retro_backend"] == "retrosim":
            pass
        return v
