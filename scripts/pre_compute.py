import argparse
import sys
import time
from configs import db_config
from multiprocessing import Pool
from pymongo import errors, MongoClient
from rdchiral.template_extractor import extract_from_reaction
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Any


def parse_args():
    parser = argparse.ArgumentParser("precompute")
    parser.add_argument("--mode",
                        help="pre-computation mode",
                        choices=["reactions"], type=str, default="reactions")

    return parser.parse_args()


def _get_products_and_fps(rxn: dict) -> dict[str, Any] | None:
    if "reaction_smiles" not in rxn:
        return None

    reaction_smiles = rxn["reaction_smiles"].split()[0]
    reactant_smiles, _, product_smiles = reaction_smiles.split(">")
    product_mol = Chem.MolFromSmiles(product_smiles)
    if product_mol is None:
        return None
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    if reactant_mol is None:
        return None

    if "reaction_smarts" in rxn:
        # user-provided
        reaction_smarts = rxn["reaction_smarts"]
    else:
        try:
            template = extract_from_reaction({
                "_id": None,
                "reactants": reactant_smiles,
                "products": product_smiles
            })
        except Exception as e:
            return None

        reaction_smarts = template.get("reaction_smarts")
        if not reaction_smarts:
            return None

    mfp_bits = list(AllChem.GetMorganFingerprintAsBitVect(
        product_mol, radius=2, nBits=2048).GetOnBits())
    pfp_bits = list(Chem.rdmolops.PatternFingerprint(
        product_mol).GetOnBits())

    mol = {
        "_id": rxn["_id"],
        "template_set": rxn["template_set"],
        "product_smiles": product_smiles,
        "reaction_smarts": reaction_smarts,
        "mfp_bits": mfp_bits,
        "mfp_count": len(mfp_bits),
        "pfp_bits": pfp_bits,
        "pfp_count": len(pfp_bits)
    }

    return mol


def precompute_reactions() -> None:
    client = MongoClient(serverSelectionTimeoutMS=1000, **db_config.MONGO)
    database = "askcos"
    collection = "reactions"
    mol_collection = "products_in_reactions"
    count_collection = "fp_counts_in_reactions"

    try:
        client.server_info()
    except errors.ServerSelectionTimeoutError:
        raise ValueError("Cannot connect to mongodb for reactions")
    else:
        db = client[database]
        db.drop_collection(mol_collection)
        db.drop_collection(count_collection)

        collection = db[collection]
        mol_collection = db[mol_collection]
        count_collection = db[count_collection]

    mfp_counts = {}
    p = Pool()
    start = time.time()
    for i, mol in enumerate(p.imap(_get_products_and_fps, collection.find())):
        if i % 100000 == 0:
            print(f"Processed {i} reactions in {time.time() - start: .2f} seconds")
            sys.stdout.flush()

        if mol is None:
            continue

        template_set = mol["template_set"]
        if template_set not in mfp_counts:
            mfp_counts[template_set] = {}
        mfp_count = mfp_counts[template_set]

        for bit in mol["mfp_bits"]:
            mfp_count[bit] = mfp_counts.get(bit, 0) + 1

        mol_collection.insert_one(mol)

    p.close()
    p.join()

    for template_set, mfp_count in mfp_counts.items():
        for k, v in mfp_count.items():
            count_collection.insert_one({
                "_id": k,
                "count": v,
                "template_set": template_set
            })

    mol_collection.create_index("mfp_bits")
    mol_collection.create_index("mfp_count")
    mol_collection.create_index("pfp_bits")
    mol_collection.create_index("pfp_count")


def main():
    args = parse_args()
    if args.mode == "reactions":
        precompute_reactions()
    else:
        raise NotImplementedError


if __name__ == "__main__":
    main()
