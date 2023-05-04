import importlib
import os
from wrappers import WRAPPER_CLASSES

_wrapper_registry = None


def get_wrapper_registry():
    """Get global wrapper registry."""
    global _wrapper_registry
    if _wrapper_registry is None:
        default_path = "configs.model_config_full"
        config_path = os.environ.get("ASKCOS_CONFIG_PATH", default_path)
        _wrapper_registry = WrapperRegistry(config_path=config_path)
        print(f"Loaded model configuration file from {config_path}")

    return _wrapper_registry


class WrapperRegistry:
    def __init__(self, config_path: str):
        config_module = importlib.import_module(config_path)
        model_config = config_module.model_config

        self._wrappers = {}
        for model_name, to_start in model_config["models_to_start"].items():
            if to_start:
                wrapper_class = WRAPPER_CLASSES[model_name]
                wrapper_config = model_config[model_name]
                self._wrappers[model_name] = wrapper_class(wrapper_config)

    def get_backend_status(self):
        return "pass"

    def get_wrapper(self, model_name: str):
        return self._wrappers[model_name]

    def __iter__(self):
        return iter(self._wrappers.values())
