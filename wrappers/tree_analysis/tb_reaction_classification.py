import networkx as nx
import traceback
from wrappers.reaction_classification.get_top_class_batch import GetTopClassBatchInput
from wrappers.registry import get_wrapper_registry
from wrappers.tree_analysis.tree_analysis_utils import (
    NIL_UUID,
    NODE_LINK_ATTRS,
    nx_paths_to_json,
    tree_data_to_graph
)


def _tb_reaction_classification(tb_result: dict) -> tuple[dict, dict]:
    """
    Run a reaction classification prediction for a saved tree builder result.

    Returns:
        dict, dict: result document and info dict with success and error fields
    """
    output = {
        "success": True,
        "error": None,
    }
    json_format = "nodelink"
    try:
        # TODO: Consolidate tree processing with duplicate code in results API
        graph = nx.node_link_graph(tb_result["graph"])  # Full reaction network
        trees = tb_result["paths"]
        if trees:
            if "nodes" in trees[0]:
                trees = [
                    nx.node_link_graph(tree, attrs=NODE_LINK_ATTRS) for tree in trees
                ]
            else:
                json_format = "treedata"
                trees = [tree_data_to_graph(tree) for tree in trees]
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "Unable to load requested result."
        print("Reaction classification failed for tree builder result:", str(e))

        return tb_result, output

    try:
        reactions = [v for v, d in graph.nodes(data=True) if d["type"] == "reaction"]
        reaction_classifier = get_wrapper_registry().get_wrapper(
            module="get_top_class_batch"
        )
        wrapper_input = GetTopClassBatchInput(smiles=reactions)
        rxn_classes = reaction_classifier.call_sync(wrapper_input).result
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "Reaction classification prediction failed."
        print("Reaction classification failed for tree builder result:", str(e))

        return tb_result, output

    try:
        # Update reaction nodes in full network
        for i, rxn in enumerate(reactions):
            rxn_data = graph.nodes[rxn]
            rxn_data["class_num"], rxn_data["class_name"] = rxn_classes[i]
        # Update reaction nodes in every tree
        for tree in trees:
            for v, d in tree.nodes(data=True):
                if d["type"] == "reaction":
                    rxn_data = graph.nodes[d["smiles"]]
                    d["class_num"], d["class_name"] = (
                        rxn_data["class_num"],
                        rxn_data["class_name"],
                    )
        # Save updated result to database
        tb_result["graph"] = nx.node_link_data(graph)
        tb_result["paths"] = nx_paths_to_json(trees, NIL_UUID, json_format=json_format)
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "Reaction classification result processing failed."
        print("Reaction classification failed for tree builder result:", str(e))

        return tb_result, output

    return tb_result, output
