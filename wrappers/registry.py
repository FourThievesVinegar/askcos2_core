import importlib
import os
import requests
from pydantic import BaseModel
from wrappers import WRAPPER_CLASSES

_wrapper_registry = None


class BackendStatus(BaseModel):
    to_start: bool
    started: bool
    backend_url: str | None


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
        self.model_config = model_config

        self._wrappers = {}
        for model_name, to_start in model_config["models_to_start"].items():
            if to_start:
                wrapper_class = WRAPPER_CLASSES[model_name]
                wrapper_config = model_config[model_name]
                self._wrappers[model_name] = wrapper_class(wrapper_config)

    def get_backend_status(self) -> dict[str, BackendStatus]:
        status = {}
        for model_name, to_start in self.model_config["models_to_start"].items():
            backend_status = {"to_start": to_start}
            wrapper = self.get_wrapper(model_name=model_name)
            if wrapper is None:
                backend_status["started"] = False
                backend_status["backend_url"] = None
            elif requests.get(wrapper.prediction_url).status_code == 404:
                backend_status["started"] = False
                backend_status["backend_url"] = None
            else:
                backend_status["started"] = True
                backend_status["backend_url"] = wrapper.prediction_url
            status[model_name] = BackendStatus(**backend_status)

        return status

    def get_wrapper(self, model_name: str):
        return self._wrappers.get(model_name, None)

    def __iter__(self):
        return iter(self._wrappers.values())
