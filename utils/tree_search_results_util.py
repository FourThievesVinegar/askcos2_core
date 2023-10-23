import networkx as nx
from rdkit import Chem
from rdkit.Chem import Descriptors


NODE_LINK_ATTRS = {
    "source": "from",
    "target": "to",
    "name": "id",
    "key": "key",
    "link": "edges",
}


def molecular_weight(smiles):
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
            return None
        mol.UpdatePropertyCache(strict=False)
        molwt = Descriptors.ExactMolWt(mol)

        return molwt


def tree_data_to_graph(tree):
    """
    Convert a single tree from tree data format to networkx graph.

    Args:
        tree (dict): tree in networkx tree data json format

    Returns:
        tree (nx.DiGraph): tree as networkx graph
    """
    # Extract attributes if they exist
    attributes = tree.pop("attributes", {})

    # Create graph from json
    tree = nx.tree_graph(tree)

    # Assign attributes to graph
    tree.graph.update(attributes)

    for node, node_data in tree.nodes.items():
        if node_data.get("is_chemical"):
            node_data["type"] = "chemical"
            del node_data["is_chemical"]
            if "terminal" not in node_data:
                # Set terminal attribute based on number of successors
                node_data["terminal"] = tree.out_degree(node) == 0
        elif node_data.get("is_reaction"):
            node_data["type"] = "reaction"
            del node_data["is_reaction"]

    return tree


def json_to_nx_paths(paths):
    """
    Convert list of paths from json to networkx graphs.
    """
    if paths:
        if "nodes" in paths[0]:
            # Nodelink format
            paths = [nx.node_link_graph(path, attrs=NODE_LINK_ATTRS) for path in paths]
        else:
            # Treedata format
            paths = [tree_data_to_graph(path) for path in paths]
    return paths


def graph_to_json(graph, strip=False):
    """
    Convert graph to json in nodelink format.
    """
    graph_json = nx.node_link_data(graph, attrs=NODE_LINK_ATTRS)
    if strip:
        keys_to_keep = {"id", "smiles", "type"}

        graph_json["nodes"] = [
            {key: value for key, value in node.items() if key in keys_to_keep}
            for node in graph_json["nodes"]
        ]
    return graph_json


def augment_v2_graph(graph):
    """
    Add missing metadata from v2 graph.

    If ``extra=True``, computes additional metadata which is
    computed at runtime by the Python tree builder, but not by
    the C++ tree builder.

    Args:
        graph (nx.Graph): networkx graph to update metadata for
        extra (optional, bool): compute extra metadata
    """
    for node, node_data in graph.nodes(data=True):
        if node_data["type"] == "chemical":
            # Calculate molecular weight
            node_data["molwt"] = molecular_weight(node)


def calculate_path_metadata(paths, graph):
    """
    Calculate pathway-level metadata to support sorting.

    Args:
        paths (list): list of nx.Graph objects corresponding to pathways
        graph (nx.Graph): data graph containing full metadata for nodes
    """
    for path in paths:
        root_id = [v for v, d in path.in_degree() if d == 0][0]
        first_reaction_id = next(path.successors(root_id))
        first_reaction_smiles = path.nodes[first_reaction_id]["smiles"]
        path.graph["first_step_score"] = graph.nodes[first_reaction_smiles][
            "rxn_score_from_model"
        ]
        path.graph["first_step_plausibility"] = graph.nodes[first_reaction_smiles][
            "plausibility"
        ]

        all_reaction_ids = (
            n for n, d in path.nodes(data=True) if d["type"] == "reaction"
        )
        all_reaction_smiles = (path.nodes[n]["smiles"] for n in all_reaction_ids)
        all_reactions = [graph.nodes[smi] for smi in all_reaction_smiles]
        num_reactions = len(all_reactions)
        path.graph["num_reactions"] = num_reactions
        path.graph["avg_score"] = (
            sum(r["rxn_score_from_model"] for r in all_reactions) / num_reactions
        )
        path.graph["avg_plausibility"] = (
            sum(r["plausibility"] for r in all_reactions) / num_reactions
        )
        path.graph["min_score"] = min(r["rxn_score_from_model"] for r in all_reactions)
        path.graph["min_plausibility"] = min(r["plausibility"] for r in all_reactions)

        if "depth" not in path.graph:
            path.graph["depth"] = [
                path.nodes[v]["type"] for v in nx.dag_longest_path(path)
            ].count("reaction")

        precursor_ids = (v for v, d in path.out_degree() if d == 0)
        precursor_smiles = [path.nodes[n]["smiles"] for n in precursor_ids]
        if "precursor_cost" not in path.graph:
            precursor_prices = [
                graph.nodes[smi]["purchase_price"] for smi in precursor_smiles
            ]
            path.graph["precursor_cost"] = (
                sum(precursor_prices) if all(precursor_prices) else None
            )

        if "atom_economy" not in path.graph:
            target = path.nodes[root_id]["smiles"]
            target_molwt = graph.nodes[target]["molwt"]
            precursor_molwt = sum(graph.nodes[smi]["molwt"]
                                  for smi in precursor_smiles)
            path.graph["atom_economy"] = target_molwt / precursor_molwt


def standardize_result_v2(result: dict) -> tuple[dict, list]:
    """
    Takes a v2 tree builder result and convert to standard format for frontend.
    """
    data_graph = result["graph"]
    paths = result["paths"]

    if paths:
        paths = json_to_nx_paths(paths)

    graph = nx.node_link_graph(data_graph)

    # Updates graph in place
    augment_v2_graph(graph)

    # Calculate some pathway metadata
    calculate_path_metadata(paths, graph)

    data_graph = nx.node_link_data(graph)
    paths = [graph_to_json(path, strip=True) for path in paths]

    return data_graph, paths


def standardize_result(result_doc: dict) -> dict:
    standardized_graph, standardized_paths = standardize_result_v2(result_doc)
    result_doc.update({
        "graph": standardized_graph,
        "paths": standardized_paths
    })

    return result_doc
