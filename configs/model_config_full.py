model_config = {
    "models_to_start": {
        "atom_map_indigo": False,
        "atom_map_rxnmapper": True
    },

    "atom_map_indigo": {
        "custom_prediction_url": "http://molgpu01.mit.edu:9918/rxnmapper",
        "timeout": 3,
        "model": {},
        "deployment": {
            "remote_image": "",
            "ports_to_expose": [],
            "start_command": ""
        },
        "celery": {
            "n_workers": 0,
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "atom_map_rxnmapper": {
        "custom_prediction_url": "http://molgpu01.mit.edu:9918/rxnmapper",
        "timeout": 3,
        "model": {},
        "deployment": {
            "remote_image": "",
            "ports_to_expose": [],
            "start_command": ""
        },
        "celery": {
            "n_workers": 0,
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    }
}
