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
        mol_collection = "products_in_reactions"
        count_collection = "counts_in_reactions"

        try:
            self.client.server_info()
        except errors.ServerSelectionTimeoutError:
            raise ValueError("Cannot connect to mongodb for reactions")
        else:
            self.db = self.client[database]
            self.collection = self.db[collection]
            self.molecules = self.db[mol_collection]
            self.counts = self.db[count_collection]

            col_present = set(self.db.list_collection_names())
            col_needed = {"reactions", "products_in_reactions", "counts_in_reactions"}

            if util_config["force_recompute_mols"]:
                self.db.drop_collection(mol_collection)
                self.db.drop_collection(count_collection)

                self.get_products_from_reactions()
                self.precompute_mols()
            elif set(col_needed) & set(col_present) != set(col_needed):
                self.db.drop_collection(mol_collection)
                self.db.drop_collection(count_collection)

                self.get_products_from_reactions()
                self.precompute_mols()

    def lookup_similar_smiles(
        self,
        smiles: str,
        sim_threshold: float,
        reaction_set: str = "USPTO_full",
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
                mol=query_mol,
                mol_collection=self.molecules,
                count_collection=self.counts,
                threshold=sim_threshold,
                reaction_set=reaction_set
            )
        elif method == "naive":
            results = sim_search(
                mol=query_mol,
                mol_collection=self.molecules,
                count_collection=self.counts,
                threshold=sim_threshold,
                reaction_set=reaction_set
            )
        elif method == "fast":
            raise NotImplementedError
        else:
            raise ValueError(f"Similarity search method '{method}' not implemented")

        output = [{'smiles': i['smiles'], 'tanimoto': i['tanimoto'], 'id': i['_id']} for i in results]

        return output

    def search_reaction_id(
        self,
        id: str,
        reaction_set: str = "USPTO_full"
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

    def get_products_from_reactions(self) -> None:
        """
        Extracts the products from each reaction in the reactions collection
        and stores them in the molecules collection to use for similarity search.
        Assumes that the database and the reactions collection have already been
        initialized.
        """
        print("Extracting products from reactions..")
        for rxn in self.collection.find():
            if "products" in rxn:
                product_smiles = rxn["products"]
            elif "reaction_smiles" in rxn:
                product_smiles = rxn["reaction_smiles"].split(">")[-1]
            else:
                continue

            self.molecules.insert_one({
                "_id": rxn["_id"],
                "smiles": product_smiles,
                "template_set": rxn["template_set"]
            })
        print("Extraction donn.")

    def precompute_mols(self) -> None:
        mfp_counts = {}
        for mol in self.molecules.find():
            smiles = mol["smiles"]
            rdmol = Chem.MolFromSmiles(smiles)
            mfp = list(AllChem.GetMorganFingerprintAsBitVect(rdmol, radius=2, nBits=2048).GetOnBits())
            pfp = list(Chem.rdmolops.PatternFingerprint(rdmol).GetOnBits())
            # smiles = Chem.MolToSmiles(rdmol)

            for bit in mfp:
                mfp_counts[bit] = mfp_counts.get(bit, 0) + 1
            
            # update the molecule document with the fingerprints and their counts
            self.molecules.update_one(
                filter={
                    "_id": mol["_id"],
                    "template_set": mol["template_set"]
                },
                update={
                    "$set": {
                        "mfp": {"bits": mfp, "count": len(mfp)},
                        "pfp": {"bits": pfp, "count": len(pfp)}
                    }
                }
            )

        for k, v in mfp_counts.items():
            self.counts.insert_one({
                "_id": k,
                "count": v,
                "template_set": mol["template_set"]
            })
        
        self.molecules.create_index("mfp.bits")
        self.molecules.create_index("mfp.count")
        self.molecules.create_index("pfp.bits")
        self.molecules.create_index("pfp.count")
