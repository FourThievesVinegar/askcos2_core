import json
from fastapi import HTTPException, Request, Response
from pydantic import BaseModel
from typing import Any, Literal
from utils import register_util
from utils.draw_impl import (
    molecule_smiles_to_image,
    reaction_smiles_to_image,
    template_smarts_to_image
)
from utils.image_annotation import generate_annotated_image


class DrawerInput(BaseModel):
    smiles: str
    input_type: Literal["chemical", "reaction", "template"] | None = None
    svg: bool | None = None
    transparent: bool | None = None
    draw_map: bool | None = None
    highlight: bool | None = None
    reacting_atoms: list | None = None
    reference: str | None = None
    align: bool | None = None
    annotate: bool | None = None
    ppg: float | None = None
    as_reactant: int | None = None
    as_product: int | None = None
    size: float | None = None


@register_util(name="draw")
class Drawer:
    """
    Util class for drawing molecules, reactions, and templates.
    Alias for legacy DrawerAPIView

    Notes:

    - Both GET and POST requests are possible. GET requests may be easier
      for simple linking, while POST requests are better for complex data.
    - If `input_type` is not specified, will attempt to determine type.
      Specifying `input_type` can provide faster results.

    Method: GET, POST

    Parameters:

    - `smiles` (str): input SMILES (or SMARTS) string
    - `input_type` (str, optional): one of 'chemical', 'reaction', or 'template'
    - `svg` (bool, optional): if True, return svg, otherwise return png
    - `transparent` (bool, optional): whether background should be transparent
    - `draw_map` (bool, optional): whether atom mapping should be drawn
    - `highlight` (bool, optional): whether to highlight mapped atoms
    - `reacting_atoms` (list, optional): list of atom scores to highlight and label
            (chemical only)
    - `reference` (str, optional): reference SMILES for aligning a molecule
    - `align` (bool, optional): whether to align reaction reactants to product
    - `annotate` (bool, optional): whether to add an annotation label to image
            (chemical only)
    - `ppg` (float, optional): price of molecule for annotation
    - `as_reactant` (int, optional): molecule popularity in training data as a
            reactant for annotation
    - `as_product` (int, optional): molecule popularity in training data as a
            product for annotation
    - `size` (float, optional): size of drawing

    Returns: SVG or PNG image of input SMILES
    """
    prefixes = ["draw"]
    methods_to_bind: dict[str, list[str]] = {
        "__call__": ["GET", "POST"]
    }

    def __init__(self, util_config: dict[str, Any]):
        pass

    def __call__(self, request: Request) -> Response:
        if request.method == "GET":
            data = dict(request.query_params)
        elif request.method == "POST":
            data = request.json()
        else:
            raise HTTPException(405, f"Unsupported method: {request.method}")

        DrawerInput(**data)         # merely validate the data

        return draw(data)


def draw(data: dict) -> Response:
    """
    Return HttpResponse with PNG of requested structure.
    """
    input_type = data.get("input_type")

    if input_type == "chemical" or ">" not in data.get("smiles"):
        methods = [draw_chemical]
    elif input_type == "reaction":
        methods = [draw_reaction]
    elif input_type == "template":
        methods = [draw_template]
    elif ">" in data.get("smiles"):
        methods = [draw_reaction, draw_template]
    else:
        methods = [draw_chemical, draw_reaction, draw_template]

    for method in methods:
        try:
            response = method(data)
        except Exception:
            response = method(data)
        else:
            return response
    else:
        d = {
            "error": "Could not draw requested structure.",
            "request": data
        }

        return Response(
            content=json.dumps(d),
            status_code=400,
            media_type="application/json"
        )


def draw_chemical(data: dict) -> Response:
    """
    Returns HttpResponse containing PNG of chemical SMILES.
    """
    smiles = data.get("smiles")
    svg = data.get("svg")
    transparent = data.get("transparent")
    highlight = data.get("highlight")
    reacting_atoms = data.get("reacting_atoms", [])
    reference = data.get("reference")
    annotate = data.get("annotate")

    media_type = "image/svg+xml" if svg else "image/png"
    if annotate:
        content = generate_annotated_image(
            smiles,
            ppg=data.get("ppg"),
            as_reactant=data.get("as_reactant"),
            as_product=data.get("as_product"),
            size=data.get("size"),
            svg=svg,
            transparent=transparent,
            highlight=highlight,
            reacting_atoms=reacting_atoms,
            reference=reference
        )
    else:
        content = molecule_smiles_to_image(
            smiles,
            svg=svg,
            transparent=transparent,
            highlight=highlight,
            reacting_atoms=reacting_atoms,
            reference=reference
        )
    response = Response(content=content, media_type=media_type)

    return response


def draw_template(data: dict) -> Response:
    """
    Returns Response containing PNG of reaction SMARTS.
    """
    smiles = data.get("smiles")
    svg = data.get("svg")
    transparent = data.get("transparent")

    media_type = "image/svg+xml" if svg else "image/png"
    content = template_smarts_to_image(smiles, svg=svg, transparent=transparent)
    response = Response(content=content, media_type=media_type)

    return response


def draw_reaction(data: dict) -> Response:
    """
    Returns Response containing PNG of reaction SMILES.
    """
    smiles = data.get("smiles")
    svg = data.get("svg")
    transparent = data.get("transparent")
    highlight = data.get("highlight")
    clear_map = not data.get("draw_map")
    align = data.get("align")

    media_type = "image/svg+xml" if svg else "image/png"
    content = reaction_smiles_to_image(
        smiles,
        svg=svg,
        transparent=transparent,
        highlight=highlight,
        clear_map=clear_map,
        align=align
    )

    response = Response(content=content, media_type=media_type)

    return response