module_config = {
    "modules_to_start": {
        "atom_map_indigo": False,
        "atom_map_rxnmapper": False,
        "forward_augmented_transformer": True,
        "general_selectivity": False
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

    "forward_augmented_transformer": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/forward_predictor/augmented_transformer.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9510, 9511, 9512],
            "default_prediction_url": "http://0.0.0.0:9510/predictions",
            "custom_prediction_url": "",
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
    }
}
