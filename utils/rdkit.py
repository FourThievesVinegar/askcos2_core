import json
import traceback as tb
from fastapi import Response
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import Descriptors, rdDepictor
from typing import Any
from utils import register_util
from utils.draw_impl import align_molecule


class SmilesInput(BaseModel):
    smiles: str
    isomericSmiles: bool = True
    reference: str = None


class MolfileInput(BaseModel):
    molfile: str
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
        "validate": ["POST"],
        "from_molfile": ["POST"],
        "to_molfile": ["POST"]
    }

    def __init__(self, util_config: dict[str, Any] = None):
        pass

    @staticmethod
    def canonicalize(input: SmilesInput) -> Response:
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
    def validate(input: SmilesInput) -> Response:
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

    @staticmethod
    def from_molfile(input: MolfileInput):
        """
        Convert the provided Molfile to a SMILES string.

        Method: POST

        Parameters:

        - `molfile` (str): Molfile input
        - `isomericSmiles` (bool, optional): whether to generate isomeric SMILES

        Returns:

        - `smiles` (str): canonical SMILES
        """
        molfile = input.molfile
        isomericSmiles = input.isomericSmiles

        resp = {}

        mol = Chem.MolFromMolBlock(molfile)
        if not mol:
            resp["error"] = "Cannot parse sdf molfile with rdkit."

            return Response(
                content=json.dumps(resp),
                status_code=500,
                media_type="application/json"
            )

        try:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
        except Exception:
            resp["error"] = "Cannot parse sdf molfile with rdkit."

            return Response(
                content=json.dumps(resp),
                status_code=500,
                media_type="application/json"
            )

        resp["smiles"] = smiles

        return Response(
            content=json.dumps(resp),
            status_code=200,
            media_type="application/json"
        )

    @staticmethod
    def to_molfile(input: SmilesInput):
        """
        Convert the provided SMILES string to a Molfile.

        Method: POST

        Parameters:

        - `smiles` (str): SMILES input
        - `reference` (str, optional): SMILES of reference molecule for alignment

        Returns:

        - `molfile` (str): Molfile output
        """

        smiles = input.smiles
        reference = input.reference

        resp = {}

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            resp["error"] = "Cannot parse smiles with rdkit."

            return Response(
                content=json.dumps(resp),
                status_code=500,
                media_type="application/json"
            )

        if reference:
            ref = Chem.MolFromSmiles(reference)
            align_molecule(mol, ref)
        else:
            rdDepictor.Compute2DCoords(mol)

        try:
            molfile = Chem.MolToMolBlock(mol)
        except Exception:
            resp["error"] = "Cannot parse smiles with rdkit."

            return Response(
                content=json.dumps(resp),
                status_code=500,
                media_type="application/json"
            )

        resp["molfile"] = molfile

        return Response(
            content=json.dumps(resp),
            status_code=200,
            media_type="application/json"
        )
