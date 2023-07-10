module_config = {
    "modules_to_start": {
        "atom_map_indigo": False,
        "atom_map_rxnmapper": False,
        "fast_filter": True,
        "forward_augmented_transformer": True,
        "forward_graph2smiles": True,
        "general_selectivity": False,
        "site_selectivity": True,
        "reaction_classification": True
    },

    "global": {
        "container_runtime": "docker",
        "image_policy": "build_all",
        "enable_gpu": True
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
    },

    "fast_filter": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/fast_filter.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9611],
            "default_prediction_url": "http://0.0.0.0:9611",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "default",
                "with_threshold"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "forward_augmented_transformer": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/forward_predictor/augmented_transformer.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9510, 9511, 9512],
            "default_prediction_url": "http://0.0.0.0:9510/predictions",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "USPTO_480k_mix"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "forward_graph2smiles": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/forward_predictor/graph2smiles.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9520, 9521, 9522],
            "default_prediction_url": "http://0.0.0.0:9520/predictions",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "USPTO_480k_mix"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "general_selectivity": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/general_selectivity.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9641]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "site_selectivity": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/site_selectivity.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9601],
            "default_prediction_url": "http://0.0.0.0:9601/site_selectivity",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "default"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "reaction_classification": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/reaction_classification.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9621],
            "default_prediction_url": "http://0.0.0.0:9621/reaction_class",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "default"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },
}
