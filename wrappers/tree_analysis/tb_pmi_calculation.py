import networkx as nx
import traceback
from wrappers.pmi_calculator.default import PmiCalculatorInput
from wrappers.registry import get_wrapper_registry
from wrappers.tree_analysis.tree_analysis_utils import (
    clean_json,
    NIL_UUID,
    NODE_LINK_ATTRS
)


def _tb_pmi_calculation(tb_result: dict, index: int) -> tuple[dict, dict]:
    if index == -1:
        paths = tb_result["paths"]
    else:
        paths = [tb_result["paths"][index]]

    output = {
        "success": True,
        "error": None,
    }

    try:
        # TODO: Consolidate tree processing with duplicate code in results API
        trees = paths
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
            output["error"] = \
                "Requested result does not have any paths for PMI calculation."

            return tb_result, output
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "Unable to load requested result."
        print("PMI calculation failed for tree builder result:", str(e))

        return tb_result, output

    try:
        # Run PMI calculator using JSON paths
        pmi_calculator = get_wrapper_registry().get_wrapper(module="pmi_calculator")
        wrapper_input = PmiCalculatorInput(trees=json_paths)
        pmis = pmi_calculator.call_sync(wrapper_input).result
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "PMI calculation failed."
        print("PMI calculation failed for tree builder result:", str(e))

        return tb_result, output

    try:
        # Update original JSON paths directly to keep same format
        for i, tree in enumerate(trees):
            # Key for graph attributes based on JSON format
            key = "graph" if "nodes" in tree else "attributes"
            tree.setdefault(key, {})["pmi"] = pmis[i]
    except Exception as e:
        traceback.print_tb(e.__traceback__)
        output["success"] = False
        output["error"] = "PMI calculation result processing failed."
        print("PMI calculation failed for tree builder result:", str(e))

        return tb_result, output

    return tb_result, output
