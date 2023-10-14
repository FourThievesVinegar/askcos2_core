import json
import traceback as tb
from fastapi import Response
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import Any
from utils import register_util


class RDKitInput(BaseModel):
    smiles: str
    isomericSmiles: bool = True


def molecular_weight(smiles: str) -> float:
    """
    Calculate exact molecular weight for the given SMILES string

    Args:
        smiles: SMILES string for which to calculate molecular weight

    Returns:
         float: exact molecular weight
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        molwt = Descriptors.ExactMolWt(mol)

        return molwt
    except Exception:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return 9999.0
        mol.UpdatePropertyCache(strict=False)
        molwt = Descriptors.ExactMolWt(mol)

        return molwt


def _canonicalize(_smi: str, isomericSmiles: bool) -> str:
    if _smi:
        _mol = Chem.MolFromSmiles(_smi)
        if not _mol:
            raise ValueError("Cannot parse smiles with rdkit.")
        _smi = Chem.MolToSmiles(_mol, isomericSmiles=isomericSmiles)
        if not _smi:
            raise ValueError("Cannot canonicalize smiles with rdkit.")
    return _smi


@register_util(name="rdkit")
class RDKitUtil:
    """Util class for RDKit"""
    prefixes = ["rdkit"]
    methods_to_bind: dict[str, list[str]] = {
        "canonicalize": ["POST"],
        "validate": ["POST"]
    }

    def __init__(self, util_config: dict[str, Any] = None):
        pass

    @staticmethod
    def canonicalize(input: RDKitInput) -> Response:
        smiles = input.smiles
        isomericSmiles = input.isomericSmiles

        resp = {}
        try:
            if ">" in smiles:
                resp["type"] = "rxn"
                resp["smiles"] = ">".join(
                    _canonicalize(part, isomericSmiles) for part in smiles.split(">")
                )
            else:
                resp["type"] = "mol"
                resp["smiles"] = _canonicalize(smiles, isomericSmiles)
        except ValueError:
            resp["error"] = f"Unable to canonicalize using RDKit, traceback: " \
                            f"{tb.format_exc()}"
            return Response(
                content=json.dumps(resp),
                status_code=500,
                media_type="application/json"
            )
        else:
            return Response(
                content=json.dumps(resp),
                status_code=200,
                media_type="application/json"
            )

    @staticmethod
    def validate(input: RDKitInput) -> Response:
        smiles = input.smiles

        mol = Chem.MolFromSmiles(smiles, sanitize=False)

        if mol is None:
            correct_syntax = False
            valid_chem_name = False
        else:
            correct_syntax = True
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                valid_chem_name = False
            else:
                valid_chem_name = True

        resp = {
            "correct_syntax": correct_syntax,
            "valid_chem_name": valid_chem_name,
        }

        return Response(
            content=json.dumps(resp),
            status_code=200,
            media_type="application/json"
        )
