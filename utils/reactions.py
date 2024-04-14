import numpy as np
import sys
import time
from configs import db_config
from pydantic import BaseModel
from pymongo import errors, MongoClient
from rdkit import Chem
from rdkit.Chem import AllChem
from schemas.base import LowerCamelAliasModel
from typing import Any
from utils import register_util
from utils.similarity_search_utils import sim_search, sim_search_aggregate

DEFAULT_MORGAN_RADIUS = 2
DEFAULT_MORGAN_LEN = 2048


class ReactionsInput(LowerCamelAliasModel):
    ids: list[int]
    template_set: str = None


class ReactionsResponse(BaseModel):
    reactions: list


def _bits_to_array(bits: list[int], fp_size: int = 2048) -> np.ndarray:
    bits_array = np.zeros(fp_size, dtype=np.bool_)
    for bit in bits:
        bits_array[bit] = True

    return bits_array


@register_util(name="reactions")
class Reactions:
    """Util class for Reactions"""
    prefixes = ["reactions"]
    methods_to_bind: dict[str, list[str]] = {
        "post": ["POST"],
        "search_reaction_id": ["POST"],
        "lookup_similar_smiles": ["POST"]
    }

    def __init__(self, util_config: dict[str, Any]):
        self.client = MongoClient(serverSelectionTimeoutMS=1000, **db_config.MONGO)
        database = "askcos"
        collection = "reactions"
        # only products with valid reaction_smarts are kept in mol_collection
        mol_collection = "products_in_reactions"
        count_collection = "fp_counts_in_reactions"

        try:
            self.client.server_info()
        except errors.ServerSelectionTimeoutError:
            raise ValueError("Cannot connect to mongodb for reactions")
        else:
            self.db = self.client[database]
            self.collection = self.db[collection]
            self.mol_collection = self.db[mol_collection]
            self.count_collection = self.db[count_collection]

            # FIXME: loading fingerprints is very expensive
            # self.mol_fps, self.mol_fp_counts, self.mol_fp_refs = (
            #     self.load_mol_fps_and_counts())

    def load_mol_fps_and_counts(self) -> tuple[
        dict[str, np.ndarray],
        dict[str, np.ndarray],
        dict[str, dict[int, str]]
    ]:
        template_sets = []
        mol_fps = {}
        mol_fp_counts = {}
        mol_fp_refs = {}

        start = time.time()
        for i, doc in enumerate(self.mol_collection.find()):
            if i % 100000 == 0:
                print(f"Processed {i} reactions in {time.time() - start: .2f} seconds")
                sys.stdout.flush()

            template_set = doc["template_set"]
            if template_set not in template_sets:
                template_sets.append(template_set)
                mol_fps[template_set] = []
                mol_fp_counts[template_set] = []
                mol_fp_refs[template_set] = {}

            current_index_in_set = len(mol_fps[template_set])
            mol_fps[template_set].append(_bits_to_array(doc["mfp_bits"]))
            mol_fp_counts[template_set].append(doc["mfp_count"])
            mol_fp_refs[template_set][current_index_in_set] = doc["_id"]

        for template_set in template_sets:
            mol_fps[template_set] = np.stack(mol_fps[template_set])
            print(f"mol_fps shape for template_set {template_set}: "
                  f"{mol_fps[template_set].shape}")
            mol_fp_counts[template_set] = np.array(
                mol_fp_counts[template_set],
                dtype=int
            )

        return mol_fps, mol_fp_counts, mol_fp_refs

    def lookup_similar_smiles(
        self,
        smiles: str,
        sim_threshold: float = 0.3,
        reaction_set: str = "USPTO_FULL",
        method: str = "accurate"
    ) -> list:
        """
        Lookup molecules in the database based on tanimoto similarity to the input
        SMILES string

        Note: assumes that a Mol Object, and Morgan Fingerprints are stored for
            each SMILES entry in the database

        Returns:
            A list of dictionary with one database entry for each molecule match including
            the tanimoto similarity to the query

        Note:
            Currently there are two options implemented lookup methods.
            The 'accurate' method is based on an aggregation pipeline in Mongo.
            The 'fast' method uses locality-sensitive hashing to greatly improve
            the lookup speed, at the cost of accuracy (especially at lower
            similarity thresholds).
        """
        query_mol = Chem.MolFromSmiles(smiles)
        if not query_mol:
            return []

        if method == "accurate":
            results = sim_search_aggregate(
                mol=query_mol,
                mol_collection=self.mol_collection,
                count_collection=self.count_collection,
                threshold=sim_threshold,
                reaction_set=reaction_set
            )
        elif method == "naive":
            results = sim_search(
                mol=query_mol,
                mol_collection=self.mol_collection,
                count_collection=self.count_collection,
                threshold=sim_threshold,
                reaction_set=reaction_set
            )
        elif method == "fast":
            raise NotImplementedError
        else:
            raise ValueError(f"Similarity search method '{method}' not implemented")

        output = [{'smiles': i['product_smiles'], 'tanimoto': i['tanimoto'], 'id': i['_id']} for i in results]

        return output

    def search_reaction_id(
        self,
        id: str,
        reaction_set: str = "USPTO_FULL"
    ) -> dict:
        """
        Lookup reaction collection using id.

        Returns:
            A dictionary containing document of reaction 
        """
        return self.collection.find_one({'_id': id, "template_set": reaction_set})

    def post(self, data: ReactionsInput) -> ReactionsResponse:
        query = {"reaction_id": {"$in": data.ids}}
        if data.template_set:
            # Processing for template subsets which use the same historian data
            query["template_set"] = data.template_set.split(":")[0]

        reactions_by_ids = list(self.collection.find(query))
        resp = ReactionsResponse(reactions=reactions_by_ids)

        return resp
