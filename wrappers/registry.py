import os

_wrapper_registry = None


def get_wrapper_registry():
    """Get global wrapper registry."""
    global _wrapper_registry
    if _wrapper_registry is None:
        default_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "model_config_full.py")
        )
        config_path = os.environ.get("ASKCOS_CONFIG_PATH", default_path)
        _wrapper_registry = WrapperRegistry(config_path=config_path)
        print(f"Loaded model configuration file from {config_path}")
    return _wrapper_registry


class WrapperRegistry:
    def __init__(self, config_path: str):
        # TODO: copy/adapt from ModelRegistry
        pass
