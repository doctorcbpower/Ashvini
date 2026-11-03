"""Ashvini: semi-analytic galaxy formation along merger trees (see README.md, MODELS.md, docs/INTERFACE_DESIGN.md)."""

# Loaded on first use, so that importing a light submodule (e.g. ashvini.reionization_physics) does not pull in the Config and run_params.yaml.
_LAZY = {"Model": "model", "PRESETS": "model", "parameters": "parameters", "TreeSample": "tree_sample"}


def __getattr__(name):
    if name in _LAZY:
        import importlib
        return getattr(importlib.import_module(f".{_LAZY[name]}", __name__), name)
    raise AttributeError(f"module 'ashvini' has no attribute {name!r}")
