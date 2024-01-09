import math
from rdkit import Chem
from rdkit.Chem import AllChem
from configs import db_config
from pydantic import BaseModel
from pymongo import errors, MongoClient
from schemas.base import LowerCamelAliasModel
from typing import Any
from utils import register_util
from utils.similarity_search_utils import sim_search, sim_search_aggregate


class ReactionsInput(LowerCamelAliasModel):
    ids: list[int]
    template_set: str = None


class ReactionsResponse(BaseModel):
    reactions: list


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

        try:
            self.client.server_info()
        except errors.ServerSelectionTimeoutError:
            raise ValueError("Cannot connect to mongodb for reactions")
        else:
            self.db = self.client[database]
            self.collection = self.client[database][collection]

            col_present = set(self.db.list_collection_names())
            col_needed = set(["molecules", "reactions", "counts"])
            if set(col_needed) & set(col_present) != set(col_needed):
                self.get_products_from_reactions()
                self.precompute_mols()
            
            self.molecules = self.db["molecules"]
            self.counts = self.db["counts"]


    def lookup_similar_smiles(
        self,
        smiles: str,
        sim_threshold: float,
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
        if method == "accurate":
            results = sim_search_aggregate(
                query_mol,
                self.molecules,
                self.counts,
                sim_threshold,
            )
        elif method == "naive":
            results = sim_search(
                query_mol,
                self.molecules,
                self.counts,
                sim_threshold,
            )
        elif method == "fast":
            raise NotImplementedError
        else:
            raise ValueError(f"Similarity search method '{method}' not implemented")

        output = [{'smiles': i['smiles'], 'tanimoto': i['tanimoto'], 'id': i['_id']} for i in results]

        return output

    def search_reaction_id(
        self,
        id: str
    ) -> dict:
        """
        Lookup reaction collection using id.

        Returns:
            A dictionary containing document of reaction 
        """
        return self.collection.find_one({'_id': id})
        # return self.molecules.find_one()
        # return self.counts.find_one()

    def post(self, data: ReactionsInput) -> ReactionsResponse:
        query = {"reaction_id": {"$in": data.ids}}
        if data.template_set:
            # Processing for template subsets which use the same historian data
            query["template_set"] = data.template_set.split(":")[0]

        reactions_by_ids = list(self.collection.find(query))
        resp = ReactionsResponse(reactions=reactions_by_ids)

        return resp

    def get_products_from_reactions(self):
        """
        Extracts the products from each reaction in the reactions collection
        and stores them in the molecules collection to use for similarity search.
        Assumes that the database and the reactions collection have already been
        initialized.
        """
        print("HELLO")
        self.molecules = self.db["molecules"]

        print("HELLO2")
        for rxn in self.collection.find():
            if rxn.get("products") is None: continue
            smiles = rxn["products"]
            id = rxn["_id"]
            self.molecules.insert_one({"_id": id, "smiles": smiles})
        print("HELLO3")

    def precompute_mols(self) -> None:
        self.counts = self.db["counts"]

        mfp_counts = {}
        for mol in self.molecules.find():
            smiles = mol["smiles"]
            rdmol = Chem.MolFromSmiles(smiles)
            mfp = list(AllChem.GetMorganFingerprintAsBitVect(rdmol, radius=2, nBits=2048).GetOnBits())
            pfp = list(Chem.rdmolops.PatternFingerprint(rdmol).GetOnBits())
            smiles = Chem.MolToSmiles(rdmol)

            for bit in mfp:
                mfp_counts[bit] = mfp_counts.get(bit, 0) + 1
            
            # update the molecule document with the fingerprints and their counts
            self.molecules.update_one({"_id": mol["_id"]}, 
                                        {"$set": {"mfp": {"bits": mfp, "count": len(mfp)},
                                                "pfp": {"bits": pfp, "count": len(pfp)}}})

        for k, v in mfp_counts.items():
            self.counts.insert_one({"_id": k, "count": v})
        
        self.molecules.create_index("mfp.bits")
        self.molecules.create_index("mfp.count")
        self.molecules.create_index("pfp.bits")
        self.molecules.create_index("pfp.count")

