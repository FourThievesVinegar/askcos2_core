model_config = {
    "models_to_start": {
        "atom_map_indigo": False,
        "atom_map_rxnmapper": True
    },

    "atom_map_indigo": {
        "timeout": 3
    },

    "atom_map_rxnmapper": {
        "custom_prediction_url": "http://molgpu01.mit.edu:9918/rxnmapper",
        "timeout": 3
    }
}
