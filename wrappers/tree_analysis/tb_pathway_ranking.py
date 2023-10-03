import networkx as nx
import traceback
from wrappers.pathway_ranker.default import PathwayRankerInput
from wrappers.registry import get_wrapper_registry
from wrappers.tree_analysis.tree_analysis_utils import (
    clean_json,
    NIL_UUID,
    NODE_LINK_ATTRS
)


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
        results = pathway_ranker.call_sync(wrapper_input).result
        results = results.dict()
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
