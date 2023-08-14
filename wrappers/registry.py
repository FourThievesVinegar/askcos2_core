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
        default_path = "configs.module_config_full"
        config_path = os.environ.get("ASKCOS_CONFIG_PATH", default_path)
        _wrapper_registry = WrapperRegistry(config_path=config_path)
        print(f"Loaded model configuration file from {config_path}")

    return _wrapper_registry


class WrapperRegistry:
    def __init__(self, config_path: str):
        module_config = importlib.import_module(config_path).module_config
        self.module_config = module_config

        self._wrappers = {}
        for module, to_start in module_config["modules_to_start"].items():
            if to_start:
                wrapper_config = module_config[module]
                if "wrapper_names" in wrapper_config:       # for multiple endpoints
                    for name in wrapper_config["wrapper_names"]:
                        wrapper_class = WRAPPER_CLASSES[name]
                        self._wrappers[name] = wrapper_class(wrapper_config)
                else:                                       # default to module name
                    wrapper_class = WRAPPER_CLASSES[module]
                    self._wrappers[module] = wrapper_class(wrapper_config)

                if module.startswith("atom_map") and \
                        "atom_map_controller" not in self._wrappers:
                    controller_class = WRAPPER_CLASSES["atom_map_controller"]
                    self._wrappers["atom_map_controller"] = controller_class()

                elif module.startswith("forward") and \
                        "forward_controller" not in self._wrappers:
                    controller_class = WRAPPER_CLASSES["forward_controller"]
                    self._wrappers["forward_controller"] = controller_class()
                #
                # if module.startswith("retro") and \
                #         "retro_controller" not in self._wrappers:
                #     controller_class = WRAPPER_CLASSES["retro_controller"]
                #     self._wrappers["retro_controller"] = controller_class

    def get_backend_status(self) -> dict[str, BackendStatus]:
        status = {}
        for module, to_start in self.module_config["modules_to_start"].items():
            backend_status = {"to_start": to_start}
            wrapper = self.get_wrapper(module=module)
            if wrapper is None:
                backend_status["started"] = False
                backend_status["backend_url"] = None
            elif requests.get(wrapper.prediction_url).status_code == 404:
                backend_status["started"] = False
                backend_status["backend_url"] = None
            else:
                backend_status["started"] = True
                backend_status["backend_url"] = wrapper.prediction_url
            status[module] = BackendStatus(**backend_status)

        return status

    def get_wrapper(self, module: str):
        return self._wrappers.get(module, None)

    def __iter__(self):
        return iter(self._wrappers.values())
