import importlib
import os
from pydantic import BaseModel

WRAPPER_CLASSES: dict[str, type] = {}


def register_wrapper(
        name: str, input_class: type[BaseModel], output_class: type[BaseModel]):
    """
    New wrapper types can be added to askcos2_core with the
    :func:`register_wrapper` function decorator.

    For example::
        @register_model(
            name='atom_map_indigo',
            input_class=AtomMapIndigoInput,
            output_class=AtomMapIndigoOutput
        )
        class AtomMapIndigoWrapper(BaseWrapper):
            (...)

    Args:
        name (str): the name of the model
        input_class (class of type[BaseModel]): the input class of the model
        output_class (class of type[BaseModel]): the output class of the model
    """

    def register_wrapper_cls(cls):
        cls.name = name
        cls.input_class = input_class
        cls.output_class = output_class

        if name not in WRAPPER_CLASSES:
            WRAPPER_CLASSES[name] = cls

        return WRAPPER_CLASSES[name]

    return register_wrapper_cls


def import_wrappers():
    wrappers_dir = os.path.dirname(__file__)
    for file in os.listdir(wrappers_dir):
        if (
            not file.startswith("_")
            and not file.startswith(".")
            and file not in ["base.py", "registry.py"]
            and (file.endswith(".py"))
        ):
            wrapper_name = file[:file.find(".py")]
            importlib.import_module(f"wrappers.{wrapper_name}")


# Automatically import any Python files in the wrappers/ directory,
# which will register all Wrapper classes with the decorator.
import_wrappers()
