import math
from rdkit.Chem import AllChem

# Code adapted from :
# https://github.com/rdkit/mongo-rdkit/blob/master/mongordkit/Search/similarity.py


# Default configurations for a variety of constants.
DEFAULT_THRESHOLD = 0.1

# Morgan fingerprints
DEFAULT_MORGAN_RADIUS = 2
DEFAULT_MORGAN_LEN = 2048

# LSH constants
DEFAULT_BIT_N = 2048


def sim_search(mol, mol_collection, count_collection=None, threshold=DEFAULT_THRESHOLD
               ) -> list:
    """
    Searches `mol_collection` for molecules with Tanimoto similarity to `mol`
    greater than or equal to `threshold`.
    Args:
        mol (rdmol Object): query molecule
        mol_collection (pymongo Collection): collection that has precomputed
            Morgan Fingerprints and bit counts.
        count_collection (pymongo Collection): popularity counts for Morgan
            fingerprint bits for all entries in `mol_collection`
        threshold: Tanimoto threshold for similarity
    Returns
        A list of MongoDB documents that fulfill `threshold`, including their
        tanimoto scores.
    """
    results = []
    qfp = list(
        AllChem.GetMorganFingerprintAsBitVect(
            mol, DEFAULT_MORGAN_RADIUS, nBits=DEFAULT_MORGAN_LEN
        ).GetOnBits()
    )

    if threshold == 0:
        for mol in mol_collection.find({}):
            tanimoto = calc_tanimoto(qfp, mol["mfp"]["bits"])
            mol["tanimoto"] = tanimoto
            results.append(mol)
        return results

    qfp_count = len(qfp)
    fp_min = int(math.ceil((threshold * qfp_count)))

    try:
        fp_max = int(qfp_count / threshold)
        print("adn here")
    except ZeroDivisionError:
        fp_max = float("inf")
    print("made it here")
    req_common_count = qfp_count - fp_min + 1
    if count_collection:
        req_common_bits = [
            count["_id"]
            for count in count_collection.find({"_id": {"$in": qfp}})
            .sort("count", 1)
            .limit(req_common_count)
        ]
        print(req_common_bits)
    else:
        req_common_bits = qfp[:req_common_count]
    for mol in mol_collection.find(
        {
            "mfp.count": {"$gte": fp_min, "$lte": fp_max},
            "mfp.bits": {"$in": req_common_bits},
        }
    ):
        tanimoto = calc_tanimoto(qfp, mol["mfp"]["bits"])
        mol["tanimoto"] = tanimoto
        if tanimoto >= threshold:
            results.append(mol)

    return results


def sim_search_aggregate(
    mol, mol_collection, count_collection=None, threshold=DEFAULT_THRESHOLD
) -> list:
    """
    Searches `mol_collection` for molecules with Tanimoto similarity to `mol`
    greater than or equal to `threshold`.
    This method uses a MongoDB aggregation pipeline.
    Args:
        mol (rdmol Object): query molecule
        mol_collection (pymongo Collection): collection that has precomputed
            Morgan Fingerprints and bit counts.
        count_collection (pymongo Collection): popularity counts for Morgan
            fingerprint bits for all entries in `mol_collection`
        threshold: Tanimoto threshold for similarity
    Returns
        A list of MongoDB documents that fulfill `threshold`, including their
        tanimoto scores.
    """
    results = []
    qfp = list(
        AllChem.GetMorganFingerprintAsBitVect(
            mol, DEFAULT_MORGAN_RADIUS, nBits=DEFAULT_MORGAN_LEN
        ).GetOnBits()
    )
    if threshold == 0:
        for mol in mol_collection.find({}):
            tanimoto = calc_tanimoto(qfp, mol["mfp"]["bits"])
            mol["tanimoto"] = tanimoto
            results.append(mol)
        return results

    qfp_count = len(qfp)
    fp_min = int(math.ceil((threshold * qfp_count)))
    try:
        fp_max = int(qfp_count / threshold)
    except ZeroDivisionError:
        fp_max = float("inf")
    req_common_count = qfp_count - fp_min + 1
    if count_collection is not None:
        req_common_bits = [
            count["_id"]
            for count in count_collection.find({"_id": {"$in": qfp}})
            .sort("count", 1)
            .limit(req_common_count)
        ]
    else:
        req_common_bits = qfp[:req_common_count]
    aggregate = [
        {
            "$match": {
                "mfp.count": {"$gte": fp_min, "$lte": fp_max},
                "mfp.bits": {"$in": req_common_bits},
            }
        },
        {
            "$set": {
                "tanimoto": {
                    "$let": {
                        "vars": {
                            "common": {
                                "$size": {"$setIntersection": ["$mfp.bits", qfp]}
                            }
                        },
                        "in": {
                            "$divide": [
                                "$$common",
                                {
                                    "$add": [
                                        qfp_count,
                                        {"$subtract": ["$mfp.count", "$$common"]},
                                    ]
                                },
                            ]
                        },
                    }
                },
            }
        },
        {"$match": {"tanimoto": {"$gte": threshold}}},
        {"$sort": {"tanimoto": -1}},
    ]
    response = mol_collection.aggregate(aggregate)

    return response


def calc_tanimoto(Na, Nb):
    """Calculates the Tanimoto similarity coefficient between two sets NA and NB."""
    Nab = len(set(Na).intersection((set(Nb))))

    return float(Nab) / (len(Na) + len(Nb) - Nab)
