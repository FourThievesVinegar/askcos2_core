import networkx as nx
import traceback as tb
from wrappers.count_analogs.default import CountAnalogsInput
from wrappers.tree_analysis.tree_analysis_utils import NODE_LINK_ATTRS
from wrappers.registry import get_wrapper_registry
from utils.registry import get_util_registry


def _tb_count_analogs_combinations(
    tree: nx.Graph,
    min_plausibility: float,
    atom_map_backend: str = "rxnmapper"
) -> int:
    """
    Run Count Combinations for calculating analogs

    Returns:
        int: result of count combinations
    """
    # FIXME together with the backend.. pricer.lookup_smarts not working
    graph_enumerator = get_wrapper_registry().get_wrapper(module="count_analogs")
    template_controller = get_util_registry().get_util(module="template")
    reaction_smiles = [
        d["smiles"] for v, d in tree.nodes(data=True) if d["type"] == "reaction"
    ]
    reaction_smarts = None
    """
    try:
        template_ids = [
            d["tforms"][0] for v, d in tree.nodes(data=True)
            if d["type"] == "reaction"
        ]
    except:
        template_ids = [
            d["id"] for v, d in tree.nodes(data=True)
            if d["type"] == "reaction" and "id" in d
        ]

    reaction_smarts = [
        template_controller.find_one_by_id(_id=t)["reaction_smarts"]
        for t in template_ids
    ]
    """
    wrapper_input = CountAnalogsInput(
        reaction_smiles=reaction_smiles,
        reaction_smarts=reaction_smarts,
        min_plausibility=min_plausibility,
        atom_map_backend=atom_map_backend
    )
    count = graph_enumerator.call_sync(wrapper_input).result

    return count


def _tb_count_analogs(
    tb_result: dict,
    index: int,
    min_plausibility: float,
    atom_map_backend: str = "rxnmapper"
) -> tuple[dict, dict]:
    """
    Run a reaction classification prediction for a saved tree builder result.

    Returns:
        dict, dict: result document and info dict with success and error fields
    """

    if index == -1:
        trees = tb_result["paths"]
    else:
        trees = [tb_result["paths"][index]]

    output = {
        "success": False,
        "error": None,
    }

    try:
        if trees:
            if "nodes" in trees[0]:
                graph_paths = [
                    nx.node_link_graph(tree, attrs=NODE_LINK_ATTRS) for tree in trees
                ]
            else:
                output["success"] = False
                output["error"] = "count_analogs only support paths in " \
                                  "nodelink format, but not treedata format"

                return tb_result, output
        else:
            output["success"] = False
            output["error"] = "Requested result does not have any paths to count."

            return tb_result, output

    except Exception as e:
        output["success"] = False
        output["error"] = "Unable to load requested result."
        print("Counting analogs failed for tree builder result:", str(e))

        return tb_result, output

    try:
        num_analogs = []
        print(f"Counting analogs for {len(graph_paths)} trees")
        for i, tree in enumerate(graph_paths):
            num_analogs.append(
                _tb_count_analogs_combinations(
                    tree=tree,
                    min_plausibility=min_plausibility,
                    atom_map_backend=atom_map_backend
                )
            )
    except Exception as e:
        output["success"] = False
        output["error"] = "Analog counting failed."
        tb.print_exc()

        return tb_result, output

    try:
        # Update target chemical nodes in every tree
        for i, tree in enumerate(trees):
            # Key for graph attributes based on JSON format
            key = "graph" if "nodes" in tree else "attributes"
            tree.setdefault(key, {})["num_analogs"] = num_analogs[i]

    except Exception as e:
        output["success"] = False
        output["error"] = "Analog counting result processing failed."
        tb.print_exc()

        return tb_result, output

    return tb_result, output
