from pydantic import BaseModel, confloat, Field
from schemas.base import LowerCamelAliasModel
from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from wrappers.legacy_solubility.utils import postprocess_solubility_results


class LegacySolubilityTask(LowerCamelAliasModel):
    solvent: str = Field(
        description="solvent SMILES"
    )
    solute: str = Field(
        description="solute SMILES"
    )
    temp: confloat(gt=0, lt=1000) = Field(
        default=298.0,
        description="temperature for solubility prediction (K)",
    )
    ref_solvent: str | None = Field(
        default=None,
        description="reference solvent SMILES"
    )
    ref_solubility: confloat(gt=-15, lt=15) | None = Field(
        default=None,
        description="solute solubility in reference solvent as logS (log10(mol/L))"
    )
    ref_temp: confloat(gt=0, lt=1000) | None = Field(
        default=None,
        description="temperature for reference solubility (K)"
    )
    hsub298: confloat(gt=0, lt=10000) | None = Field(
        default=None,
        description="sublimation enthalpy of solute at 298 K (kcal/mol)"
    )
    cp_gas_298: confloat(gt=0, lt=10000) | None = Field(
        default=None,
        description="gas phase heat capacity of solute at 298 K (cal/K/mol)"
    )
    cp_solid_298: confloat(gt=0, lt=10000) | None = Field(
        default=None,
        description="solid phase heat capacity of solute at 298 K (cal/K/mol)"
    )


class LegacySolubilityInput(LowerCamelAliasModel):
    task_list: list[LegacySolubilityTask]

    class Config:
        schema_extra = {
            "example": {
                "task_list": [
                    {
                        "solvent": "CC(=O)O",
                        "solute": "C(CCC(=O)O)CC(=O)O",
                        "temp": 298.0,
                    },
                    {
                        "solvent": "CC(=O)O",
                        "solute": "C(CCC(=O)O)CC(=O)O",
                        "temp": 500.0,
                    },
                ]
            }
        }


class LegacySolubilityOutput(BaseModel):
    Solvent: list[str]
    Solute: list[str]
    Temp: list[float]
    ref_solvent: list[str | None] = Field(alias="Ref. Solv")
    ref_solubility: list[str | float | None] = Field(alias="Ref. Solub")
    ref_temp: list[str | float | None] = Field(alias="Ref. Temp")
    hsub298: list[str | float | None] = Field(alias="Input Hsub298")
    cp_gas_298: list[str | float | None] = Field(alias="Input Cpg298")
    cp_solid_298: list[str | float | None] = Field(alias="Input Cps298")
    error_message: list[str | None] = Field(alias="Error Message")
    warning_message: list[str | None] = Field(alias="Warning Message")
    log_st_1: list[str | float | None] = Field(alias="logST (method1) [log10(mol/L)]")
    log_st_2: list[str | float | None] = Field(alias="logST (method2) [log10(mol/L)]")
    dg_solv_t: list[str | float | None] = Field(alias="dGsolvT [kcal/mol]")
    dh_solv_t: list[str | float | None] = Field(alias="dHsolvT [kcal/mol]")
    ds_solv_t: list[str | float | None] = Field(alias="dSsolvT [cal/K/mol]")
    pred_hsub298: list[str | float | None] = Field(alias="Pred. Hsub298 [kcal/mol]")
    pred_cpg298: list[str | float | None] = Field(alias="Pred. Cpg298 [cal/K/mol]")
    pred_cps298: list[str | float | None] = Field(alias="Pred. Cps298 [cal/K/mol]")
    log_s_298: list[str | float | None] = Field(alias="logS298 [log10(mol/L)]")
    uncertainty_log_s_298: list[str | float | None] = Field(
        alias="uncertainty logS298 [log10(mol/L)]")
    dh_solv_298: list[str | float | None] = Field(alias="dHsolv298 [kcal/mol]")
    uncertainty_dh_solv_298: list[str | float | None] = Field(
        alias="uncertainty dHsolv298 [kcal/mol]")
    E: list[str | float | None]
    S: list[str | float | None]
    A: list[str | float | None]
    B: list[str | float | None]
    L: list[str | float | None]
    V: list[str | float | None]


class LegacySolubilityResult(BaseModel):
    Solvent: str
    Solute: str
    Temp: float
    ref_solvent: str | None = Field(alias="Ref. Solv")
    ref_solubility: str | float | None = Field(alias="Ref. Solub")
    ref_temp: str | float | None = Field(alias="Ref. Temp")
    hsub298: str | float | None = Field(alias="Input Hsub298")
    cp_gas_298: str | float | None = Field(alias="Input Cpg298")
    cp_solid_298: str | float | None = Field(alias="Input Cps298")
    error_message: str | None = Field(alias="Error Message")
    warning_message: str | None = Field(alias="Warning Message")
    log_st_1: str | float | None = Field(alias="logST (method1) [log10(mol/L)]")
    log_st_2: str | float | None = Field(alias="logST (method2) [log10(mol/L)]")
    dg_solv_t: str | float | None = Field(alias="dGsolvT [kcal/mol]")
    dh_solv_t: str | float | None = Field(alias="dHsolvT [kcal/mol]")
    ds_solv_t: str | float | None = Field(alias="dSsolvT [cal/K/mol]")
    pred_hsub298: str | float | None = Field(alias="Pred. Hsub298 [kcal/mol]")
    pred_cpg298: str | float | None = Field(alias="Pred. Cpg298 [cal/K/mol]")
    pred_cps298: str | float | None = Field(alias="Pred. Cps298 [cal/K/mol]")
    log_s_298: str | float | None = Field(alias="logS298 [log10(mol/L)]")
    uncertainty_log_s_298: str | float | None = Field(
        alias="uncertainty logS298 [log10(mol/L)]")
    dh_solv_298: str | float | None = Field(alias="dHsolv298 [kcal/mol]")
    uncertainty_dh_solv_298: str | float | None = Field(
        alias="uncertainty dHsolv298 [kcal/mol]")
    E: str | float | None
    S: str | float | None
    A: str | float | None
    B: str | float | None
    L: str | float | None
    V: str | float | None

    st_1: str | float | None = Field(alias="ST (method1) [mg/mL]")
    st_2: str | float | None = Field(alias="ST (method2) [mg/mL]")
    s_298: str | float | None = Field(alias="S298 [mg/mL]")


class LegacySolubilityResponse(BaseModel):
    __root__: list[LegacySolubilityResult]


@register_wrapper(
    name="legacy_solubility",
    input_class=LegacySolubilityInput,
    output_class=LegacySolubilityOutput,
    response_class=LegacySolubilityResponse
)
class LegacySolubilityWrapper(BaseWrapper):
    """Wrapper class for Legacy Solubility"""
    prefixes = ["legacy/solubility/batch"]

    def call_raw(self, input: LegacySolubilityInput) -> LegacySolubilityOutput:
        input_dict = {}
        for k in LegacySolubilityTask.__fields__.keys():
            input_dict[f"{k}_list"] = [getattr(task, k) for task in input.task_list]

        response = self.session_sync.post(
            self.prediction_url,
            json=input_dict,
            timeout=self.config["deployment"]["timeout"]
        )
        output = response.json()
        output = LegacySolubilityOutput(**output)

        return output

    def call_sync(self, input: LegacySolubilityInput) -> LegacySolubilityResponse:
        """
        Endpoint for synchronous call to solubility prediction backend
        """
        output = self.call_raw(input=input)
        response = self.convert_output_to_response(output)

        return response

    async def call_async(self, input: LegacySolubilityInput, priority: int = 0) -> str:
        """
        Endpoint for asynchronous call to solubility prediction backend
        """
        return await super().call_async(input=input, priority=priority)

    async def retrieve(self, task_id: str) -> LegacySolubilityResponse | None:
        return await super().retrieve(task_id=task_id)

    @staticmethod
    def convert_output_to_response(output: LegacySolubilityOutput
                                   ) -> LegacySolubilityResponse:
        keys, value_lists = zip(*output.dict(by_alias=True).items())
        results = [dict(zip(keys, values)) for values in zip(*value_lists)]
        results = postprocess_solubility_results(results)
        response = LegacySolubilityResponse(__root__=results)

        return response
