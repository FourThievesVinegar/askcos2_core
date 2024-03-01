module_config = {
    "modules_to_start": {
        "atom_map_indigo": False,
        "atom_map_rxnmapper": True,
        "atom_map_wln": False,
        "cluster": True,
        "context_recommender": False,
        "count_analogs": True,
        "descriptors": False,
        "evaluate_reactions": False,
        "fast_filter": True,
        "forward_augmented_transformer": False,
        "forward_graph2smiles": False,
        "forward_wldn5": False,
        "general_selectivity": False,
        "impurity_predictor": False,
        "pathway_ranker": True,
        "pmi_calculator": True,
        "reaction_classification": True,
        "retro_template_relevance": True,
        "retro_augmented_transformer": False,
        "scscore": True,
        "site_selectivity": False,
        "tree_search_expand_one": True,
        "tree_search_mcts": True,
        "tree_search_retro_star": True,
        "legacy_descriptors": False,
        "legacy_solubility": False
    },

    "global": {
        "require_frontend": False,
        "container_runtime": "docker",
        "image_policy": "build_all",
        "enable_gpu": False
    },

    "atom_map_indigo": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/atom_map/indigo.git",
        "description": "Atom mapper from the indigo package",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9661],
            "default_prediction_url": "http://0.0.0.0:9661/indigo_mapper",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        }
    },

    "atom_map_rxnmapper": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/atom_map/rxnmapper.git",
        "description": "Wrapper around IBM RXNMapper",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9671],
            "default_prediction_url": "http://0.0.0.0:9671/ibm_rxnmapper",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        }
    },

    "atom_map_wln": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/atom_map/wln.git",
        "description": "Atom mapper based on the WLN model",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9651],
            "default_prediction_url": "http://0.0.0.0:9651/wln_mapper",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        }
    },

    "cluster": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/cluster.git",
        "description":
            "Clusterer using kmeans, hdbscan, or based on reaction classification",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9801],
            "default_prediction_url": "http://0.0.0.0:9801/cluster",
            "custom_prediction_url": "",
            "timeout": 30,
            "available_model_names": []
        }
    },

    "context_recommender": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/context_recommender.git",
        "description": "Context recommender trained using fingerprint or GNN",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9901],
            "default_prediction_url": "http://0.0.0.0:9901",
            "custom_prediction_url": "",
            "timeout": 20,
            "available_model_names": []
        },
        "wrapper_names": [
            "context_recommender_cleaned",
            "context_recommender_uncleaned",
            "context_recommender_fp",
            "context_recommender_graph",
            "context_recommender_pr_fp",
            "context_recommender_pr_graph"
        ]
    },

    "count_analogs": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/count_analogs.git",
        "description": "Analog counting via graph enumeration",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9911],
            "default_prediction_url": "http://0.0.0.0:9911/count_analogs",
            "custom_prediction_url": "",
            "timeout": 120,
            "available_model_names": []
        }
    },

    "descriptors": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/descriptors.git",
        "description": "Predictor for quantum descriptors",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9631],
            "default_prediction_url": "http://0.0.0.0:9631/descriptors",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        }
    },

    "evaluate_reactions": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/evaluate_reactions.git",
        "description": "Reaction evaluator for reactions with different contexts",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9721],
            "default_prediction_url": "http://0.0.0.0:9721",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        },
        "wrapper_names": [
            "evaluate_context",
            "evaluate_reaction"
        ]
    },

    "fast_filter": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/fast_filter.git",
        "description":
            "Fast binary classification model for filtering out unlikely reaction",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9611],
            "default_prediction_url": "http://0.0.0.0:9611",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        },
        "wrapper_names": [
            "fast_filter",
            "fast_filter_with_threshold",
            "fast_filter_batch"
        ]
    },

    "forward_augmented_transformer": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/forward_predictor/augmented_transformer.git",
        "description": "Forward predictor based on the Augmented Transformer model",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9510, 9511, 9512],
            "default_prediction_url": "http://0.0.0.0:9510/predictions",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": [
                "pistachio_23Q3",
                "USPTO_STEREO"
            ]
        }
    },

    "forward_graph2smiles": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/forward_predictor/graph2smiles.git",
        "description": "Forward predictor based on the Graph2SMILES model",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9520, 9521, 9522],
            "default_prediction_url": "http://0.0.0.0:9520/predictions",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": [
                "pistachio_23Q3",
                "USPTO_STEREO"
            ]
        }
    },

    "forward_wldn5": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/forward_predictor/wldn5.git",
        "description": "Forward predictor based on the WLDN5 model",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9501],
            "default_prediction_url": "http://0.0.0.0:9501/wldn5_predict",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": [
                "uspto_500k",
                "pistachio"
            ]
        }
    },

    "general_selectivity": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/general_selectivity.git",
        "description":
            "General selectivity predictor based on GNN, optionally with QM features",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9641],
            "default_prediction_url": "http://0.0.0.0:9641",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        },
        "wrapper_names": [
            "general_selectivity_gnn",
            "general_selectivity_qm_gnn",
            "general_selectivity_qm_gnn_no_reagent",
        ]
    },

    "impurity_predictor": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/impurity_predictor.git",
        "description":
            "Impurity predictor based on forward predictor(s), "
            "e.g., by exploring over-reaction",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9691],
            "default_prediction_url": "http://0.0.0.0:9691/impurity",
            "custom_prediction_url": "",
            "timeout": 180,
            "available_model_names": []
        }
    },

    "pathway_ranker": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/pathway_ranker.git",
        "description": "Reaction pathway ranking model (with Tree-LSTM)",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9681],
            "default_prediction_url": "http://0.0.0.0:9681/pathway_ranker",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        }
    },

    "pmi_calculator": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/pmi_calculator.git",
        "description": "Calculator for Process Mass Intensity",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9701],
            "default_prediction_url": "http://0.0.0.0:9701/pmi_calculator",
            "custom_prediction_url": "",
            "timeout": 120,
            "available_model_names": []
        }
    },

    "reaction_classification": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/reaction_classification.git",
        "description": "Reaction classifier trained with Pistachio dataset",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9621],
            "default_prediction_url": "http://0.0.0.0:9621",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        },
        "wrapper_names": [
            "reaction_classification",
            "get_top_class_batch"
        ]
    },

    "retro_template_relevance": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/retro/template_relevance.git",
        "description":
            "One-step retrosynthesis model with an FFN-based template classifier",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9410, 9411, 9412],
            "default_prediction_url": "http://0.0.0.0:9410/predictions",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": [
                "bkms_metabolic",
                "pistachio",
                "pistachio_ringbreaker",
                "reaxys",
                "reaxys_biocatalysis"
            ]
        }
    },

    "retro_augmented_transformer": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/retro/augmented_transformer.git",
        "description":
            "One-step retrosynthesis model using the Augmented Transformer model",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9420, 9421, 9422],
            "default_prediction_url": "http://0.0.0.0:9420/predictions",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": [
                "pistachio_23Q3",
                "USPTO_FULL"
            ]
        }
    },

    "scscore": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/scscore.git",
        "description": "SCScorer for calculating synthetic complexity",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9741],
            "default_prediction_url": "http://0.0.0.0:9741/scscore",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        },
        "wrapper_names": [
            "scscore",
            "scscore_batch"
        ]
    },

    "site_selectivity": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/site_selectivity.git",
        "description": "Site selectivity predictor",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9601],
            "default_prediction_url": "http://0.0.0.0:9601/site_selectivity",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        }
    },

    "tree_search_expand_one": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/tree_search/expand_one.git",
        "description": "The controller for one-step expansion in IPP and tree builder",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9301],
            "default_prediction_url": "http://0.0.0.0:9301/get_outcomes",
            "custom_prediction_url": "",
            "timeout": 30,
            "available_model_names": []
        }
    },

    "tree_search_mcts": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/tree_search/mcts.git",
        "description":
            "The controller for Monte Carlo Tree Search (i.e., the Tree Builder)",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9311],
            "default_prediction_url": "http://0.0.0.0:9311/get_buyable_paths",
            "custom_prediction_url": "",
            "timeout": 1200,
            "available_model_names": []
        }
    },

    "tree_search_retro_star": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/tree_search/retro_star.git",
        "description":
            "The controller for the Retro* planner",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "use_gpu": False,
            "ports_to_expose": [9321],
            "default_prediction_url": "http://0.0.0.0:9321/get_buyable_paths",
            "custom_prediction_url": "",
            "timeout": 1200,
            "available_model_names": []
        }
    },

    "legacy_descriptors": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/legacy_descriptors.git",
        "description": "Black-box quantum descriptor model from ASKCOSv1",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "image_policy": "pull",
            "use_gpu": False,
            "ports_to_expose": [9632],
            "default_prediction_url": "http://0.0.0.0:9632/predictions/descriptors",
            "custom_prediction_url": "",
            "timeout": 10,
            "available_model_names": []
        }
    },

    "legacy_solubility": {
        "repo": "git@gitlab.com:mlpds_mit/askcosv2/legacy_solubility.git",
        "description": "Black-box solubility model from ASKCOSv1",
        "deployment": {
            "deployment_config": "deployment.yaml",
            "image_policy": "pull",
            "use_gpu": False,
            "ports_to_expose": [9732],
            "default_prediction_url": "http://0.0.0.0:9732/predictions/solprop",
            "custom_prediction_url": "",
            "timeout": 60,
            "available_model_names": []
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

    "celery_task": {},

    "draw": {},

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

    "selectivity_refs": {
        "file": "utils/site_selectivity_refs.json"
    },

    "tree_search_results_controller": {
        "database": "results",
        "collection": "results"
    },

    "user_controller": {
        "engine": "db",
        "database": "askcos",
        "collection": "users"
    }
}
