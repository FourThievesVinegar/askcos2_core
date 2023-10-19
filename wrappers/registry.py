import importlib
import os
from pydantic import BaseModel
from wrappers import WRAPPER_CLASSES

_wrapper_registry = None


class BackendStatus(BaseModel):
    name: str
    description: str
    available_model_names: list[str]
    to_start: bool
    ready: bool
    backend_url: str | None


class BackendStatusResponse(BaseModel):
    modules: list[BackendStatus]


def get_wrapper_registry():
    """Get global wrapper registry."""
    global _wrapper_registry
    if _wrapper_registry is None:
        default_path = "configs.module_config_full"
        config_path = os.environ.get("ASKCOS_CONFIG_PATH", default_path)
        _wrapper_registry = WrapperRegistry(config_path=config_path)
        print(f"Loaded wrapper configuration from {config_path}")

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

                elif module.startswith("general_selectivity") and \
                        "general_selectivity_controller" not in self._wrappers:
                    controller_class = \
                        WRAPPER_CLASSES["general_selectivity_controller"]
                    self._wrappers["general_selectivity_controller"] = \
                        controller_class()

                elif module.startswith("retro") and \
                        "retro_controller" not in self._wrappers:
                    controller_class = WRAPPER_CLASSES["retro_controller"]
                    self._wrappers["retro_controller"] = controller_class()

        # tree analysis controller will be created always
        controller_class = WRAPPER_CLASSES["tree_analysis_controller"]
        self._wrappers["tree_analysis_controller"] = controller_class()

    def get_backend_status(self) -> BackendStatusResponse:
        """Endpoint for getting the status of backend services"""
        status_list = []

        for module, to_start in self.module_config["modules_to_start"].items():
            backend_url = ""
            ready = False

            if to_start:
                wrapper_config = self.module_config[module]
                if "wrapper_names" in wrapper_config:       # for multiple endpoints
                    wrapper = self.get_wrapper(
                        module=wrapper_config["wrapper_names"][0]
                    )
                else:
                    wrapper = self.get_wrapper(module=module)
                if wrapper is not None:
                    backend_url = wrapper.config["prediction_url_in_use"]
                    ready = wrapper.is_ready()

            backend_status = {
                "name": module,
                "description": self.module_config[module]["description"],
                "available_model_names":
                    self.module_config[module]["deployment"]["available_model_names"],
                "to_start": to_start,
                "ready": ready,
                "backend_url": backend_url
            }
            backend_status = BackendStatus(**backend_status)
            status_list.append(backend_status)

        resp = BackendStatusResponse(modules=status_list)

        return resp

    def get_wrapper(self, module: str):
        return self._wrappers.get(module, None)

    def __iter__(self):
        return iter(self._wrappers.values())
