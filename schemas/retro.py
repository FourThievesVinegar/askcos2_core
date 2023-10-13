from pydantic import Field, validator
from schemas.base import LowerCamelAliasModel
from typing import Any, Literal


class RetroBackendOption(LowerCamelAliasModel):
    retro_backend: Literal[
        "template_relevance",
        "augmented_transformer",
        "graph2smiles"
    ] = "template_relevance"
    retro_model_name: str = "reaxys"
    max_num_templates: int = 100
    max_cum_prob: float = 0.995
    attribute_filter: list[dict[str, Any]] = Field(default_factory=list)

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
                "pistachio_2023Q2"
            ]:
                raise ValueError(
                    f"Unsupported retro_model_name {v} for augmented_transformer")
        elif values["retro_backend"] == "graph2smiles":
            if v not in [
                "cas",
                "pistachio_2023Q2"
            ]:
                raise ValueError(
                    f"Unsupported retro_model_name {v} for graph2smiles")
        return v
