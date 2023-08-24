module_config = {
    "modules_to_start": {
        "atom_map_indigo": True,
        "atom_map_rxnmapper": True,
        "atom_map_wln": True,
        "cluster": True,
        "context_recommender": True,
        "descriptors": True,
        "evaluate_reactions": True,
        "fast_filter": True,
        "forward_augmented_transformer": True,
        "forward_graph2smiles": True,
        "forward_wldn5": True,
        "general_selectivity": True,
        "impurity_predictor": True,
        "pathway_ranker": True,
        "pmi_calculator": True,
        "reaction_classification": True,
        "retro_template_relevance": True,
        "retro_augmented_transformer": True,
        "retro_graph2smiles": True,
        "scscore": True,
        "site_selectivity": True,
        "tree_optimizer": True
    },

    "global": {
        "container_runtime": "docker",
        "image_policy": "build_all",
        "enable_gpu": True
    },

    "atom_map_indigo": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/atom_map/indigo.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9661],
            "default_prediction_url": "http://0.0.0.0:9661/indigo_mapper",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": []
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "atom_map_rxnmapper": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/atom_map/rxnmapper.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9671],
            "default_prediction_url": "http://0.0.0.0:9671/ibm_rxnmapper",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": []
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "atom_map_wln": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/atom_map/wln.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9651],
            "default_prediction_url": "http://0.0.0.0:9651/wln_mapper",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": []
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "cluster": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/cluster.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9801],
            "default_prediction_url": "http://0.0.0.0:9801/cluster",
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

    "context_recommender": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/context_recommender.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9901],
            "default_prediction_url": "http://0.0.0.0:9901",
            "custom_prediction_url": "",
            "timeout": 20,
            "available_models": [
                "default"
            ]
        },
        "wrapper_names": [
            "context_recommender_cleaned",
            "context_recommender_uncleaned",
            "context_recommender_fp",
            "context_recommender_graph",
            "context_recommender_pr_fp",
            "context_recommender_pr_graph"
        ],
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "descriptors": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/descriptors.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9631],
            "default_prediction_url": "http://0.0.0.0:9631/descriptors",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": []
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "evaluate_reactions": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/evaluate_reactions.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9721],
            "default_prediction_url": "http://0.0.0.0:9721",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_models": []
        },
        "wrapper_names": [
            "evaluate_context",
            "evaluate_reaction"
        ],
        "celery": {
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
            "available_models": []
        },
        "wrapper_names": [
            "fast_filter",
            "fast_filter_with_threshold"
        ],
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

    "forward_wldn5": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/forward_predictor/wldn5.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9501],
            "default_prediction_url": "http://0.0.0.0:9501/wldn5_predict",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "uspto_500k",
                "pistachio"
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
            "ports_to_expose": [9641],
            "default_prediction_url": "http://0.0.0.0:9641",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_models": [
                "default"
            ]
        },
        "wrapper_names": [
            "general_selectivity_gnn",
            "general_selectivity_qm_gnn",
            "general_selectivity_qm_gnn_no_reagent",
        ],
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "impurity_predictor": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/impurity_predictor.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9691],
            "default_prediction_url": "http://0.0.0.0:9691/impurity",
            "custom_prediction_url": "",
            "timeout": 60,
            "available_models": [
                "default"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "pathway_ranker": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/pathway_ranker.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9681],
            "default_prediction_url": "http://0.0.0.0:9681/pathway_ranker",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_models": [
                "default"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "pmi_calculator": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/pmi_calculator.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9701],
            "default_prediction_url": "http://0.0.0.0:9701/pmi_calculator",
            "custom_prediction_url": "",
            "timeout": 120,
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
            "default_prediction_url": "http://0.0.0.0:9621",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_models": [
                "default"
            ]
        },
        "wrapper_names": [
            "reaction_classification",
            "get_top_class_batch"
        ],
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "retro_template_relevance": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/retro/template_relevance.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9410, 9411, 9412],
            "default_prediction_url": "http://0.0.0.0:9410/predictions",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "placeholder"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "retro_augmented_transformer": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/retro/augmented_transformer.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": True,
            "ports_to_expose": [9420, 9421, 9422],
            "default_prediction_url": "http://0.0.0.0:9420/predictions",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "placeholder"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "retro_graph2smiles": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/retro/graph2smiles.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": True,
            "ports_to_expose": [9430, 9431, 9432],
            "default_prediction_url": "http://0.0.0.0:9430/predictions",
            "custom_prediction_url": "",
            "timeout": 3,
            "available_models": [
                "placeholder"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    "scscore": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/scscore.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9741],
            "default_prediction_url": "http://0.0.0.0:9741/scscore",
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

    "tree_optimizer": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/tree_optimizer.git",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9711],
            "default_prediction_url": "http://0.0.0.0:9711/tree_optimizer",
            "custom_prediction_url": "",
            "timeout": 600,
            "available_models": [
                "default"
            ]
        },
        "celery": {
            "queue_name": "generic",
            "worker_name": "generic_worker"
        }
    },

    # utils
    "banned_chemicals": {
        "engine": "db",
        "database": "askcos",
        "collection": "banned_chemicals"
    },

    "banned_reactions": {
        "engine": "db",
        "database": "askcos",
        "collection": "banned_reactions"
    },

    "historian": {
        "engine": "db",
        "database": "askcos",
        "collection": "chemicals",
        "files": []
    },

    "pricer": {
        "engine": "db",
        "database": "askcos",
        "collection": "buyables",
        "file": "",
        "precompute_mols": False
    },

    "user_controller": {
        "engine": "db",
        "database": "askcos",
        "collection": "users"
    }
}
