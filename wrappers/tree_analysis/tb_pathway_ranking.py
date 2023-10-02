import networkx as nx
import traceback
from wrappers.pathway_ranker.default import PathwayRankerInput
from wrappers.registry import get_wrapper_registry

NIL_UUID = "00000000-0000-0000-0000-000000000000"
NODE_LINK_ATTRS = {
    "source": "from",
    "target": "to",
    "name": "id",
    "key": "key",
    "link": "edges",
}
# Map from keys used by tree builder graph to name used in pathways
PATH_KEY_DICT = {
    "smiles": "smiles",
    "type": "type",
    "id": "id",
    "as_reactant": "as_reactant",
    "as_product": "as_product",
    "plausibility": "plausibility",
    "forward_score": "forward_score",   # From graph optimization
    "ppg": "ppg",           # If calling clean_json on a previously cleaned tree
    "purchase_price": "ppg",
    "properties": "properties",
    "template_score": "template_score",
    "terminal": "terminal",
    "tforms": "tforms",
    "tsources": "tsources",
    "num_examples": "num_examples",
    "necessary_reagent": "necessary_reagent",
    "precursor_smiles": "precursor_smiles",
    "rms_molwt": "rms_molwt",
    "num_rings": "num_rings",
    "scscore": "scscore",
    "rank": "rank",
    "class_num": "class_num",
    "class_name": "class_name",
}
# List of all keys to include in output JSON
OUTPUT_KEYS = {
    "chemical": [
        "smiles",
        "id",
        "as_reactant",
        "as_product",
        "ppg",
        "properties",
        "terminal",
    ],
    "reaction": [
        "smiles",
        "id",
        "plausibility",
        "forward_score",
        "template_score",
        "tforms",
        "tsources",
        "num_examples",
        "necessary_reagent",
        "precursor_smiles",
        "rms_molwt",
        "num_rings",
        "scscore",
        "rank",
        "class_num",
        "class_name",
    ],
}


def clean_json(path):
    """
    Clean up json representation of a pathway. Accepts paths from either
    tree builder version.

    Note about chemical/reaction node identification:
        * For treedata format, chemical nodes have an ``is_chemical`` attribute,
          while reaction nodes have an ``is_reaction`` attribute
        * For nodelink format, all nodes have a ``type`` attribute, whose value
          is either ``chemical`` or ``reaction``

    This distinction in the JSON schema is for historical reasons and is not
    an official aspect of the treedata/nodelink formats.
    """

    def _add_missing_fields(_node, _type):
        for _key in OUTPUT_KEYS[_type]:
            if _key not in _node:
                _node[_key] = None

    if "nodes" in path:
        # Node link format
        nodes = []
        for node in path["nodes"]:
            new_node = {
                PATH_KEY_DICT[key]: value
                for key, value in node.items()
                if key in PATH_KEY_DICT
            }
            # Set any output keys not present to None
            _add_missing_fields(new_node, node["type"])
            nodes.append(new_node)
        path["nodes"] = nodes
        output = path
    else:
        # Tree data format
        output = {}
        for key, value in path.items():
            if key == "type":
                if value == "chemical":
                    output["is_chemical"] = True
                elif value == "reaction":
                    output["is_reaction"] = True
            elif key == "children":
                output["children"] = [clean_json(c) for c in value]
            elif key in PATH_KEY_DICT:
                output[PATH_KEY_DICT[key]] = value

        # Set any output keys not present to None
        _add_missing_fields(output, path["type"])

        if "children" not in output:
            output["children"] = []

    return output


def _tb_pathway_ranking(
    tb_result: dict,
    clustering: bool,
    cluster_method: str,
    min_samples: int,
    min_cluster_size: int
) -> tuple[dict, dict]:
    """
    Run a pathway ranking prediction for a saved tree builder result.

    Args:
        tb_result (dict): tree builder saved result document

    Returns:
        dict, dict: result document and info dict with success and error fields
    """
    output = {
        "success": True,
        "error": None,
    }

    try:
        # TODO: Consolidate tree processing with duplicate code in results API
        trees = tb_result["paths"]
        if trees:
            if "nodes" in trees[0]:
                # Convert paths to tree data format for pathway ranking
                graph_paths = [
                    nx.node_link_graph(tree, attrs=NODE_LINK_ATTRS) for tree in trees
                ]
                json_paths = [
                    clean_json(nx.tree_data(tree, NIL_UUID)) for tree in graph_paths
                ]
            else:
                json_paths = trees
        else:
            output["success"] = False
            output["error"] = "Requested result does not have any paths to rank."

            return tb_result, output
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "Unable to load requested result."
        print("Pathway ranking failed for tree builder result:", str(e))

        return tb_result, output

    try:
        # Run pathway ranking and clustering using JSON paths
        pathway_ranker = get_wrapper_registry().get_wrapper(module="pathway_ranker")
        wrapper_input = PathwayRankerInput(
            trees=json_paths,
            clustering=clustering,
            cluster_method=cluster_method,
            min_samples=min_samples,
            min_cluster_size=min_cluster_size
        )
        results = pathway_ranker(wrapper_input).result
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "Pathway ranking prediction failed."
        print("Pathway ranking failed for tree builder result:", str(e))

        return tb_result, output

    try:
        # Update original JSON paths directly to keep same format
        for i, tree in enumerate(trees):
            # Key for graph attributes based on JSON format
            key = "graph" if "nodes" in tree else "attributes"
            tree.setdefault(key, {})["score"] = results["scores"][i]
            if clustering:
                tree[key]["cluster_id"] = results["clusters"][i]
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "Pathway ranking result processing failed."
        print("Pathway ranking failed for tree builder result:", str(e))

        return tb_result, output

    return tb_result, output
