import importlib
import os
from utils import UTIL_CLASSES

_util_registry = None


def get_util_registry():
    """Get global util registry."""
    global _util_registry
    if _util_registry is None:
        default_path = "configs.module_config_full"
        config_path = os.environ.get("ASKCOS_CONFIG_PATH", default_path)
        _util_registry = UtilRegistry(config_path=config_path)
        print(f"Loaded module configuration file from {config_path}")

    return _util_registry


class UtilRegistry:
    def __init__(self, config_path: str):
        module_config = importlib.import_module(config_path).module_config
        self.module_config = module_config

        self._utils = {
            util_name: util_class(util_config=module_config[util_name])
            for util_name, util_class in UTIL_CLASSES.items()
        }

    def get_util(self, module: str):
        return self._utils.get(module, None)

    def __iter__(self):
        return iter(self._utils.values())
