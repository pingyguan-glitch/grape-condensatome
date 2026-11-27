"""PlantLLPS — plant-aware LLPS propensity scoring (prototype).

Honest scope: this package ships a runnable feature pipeline, a tiny
bundled training set, a transparent composite + logistic stub model,
and a GrapeCondensatome SQLite schema. Full proteome scoring and
production weights trained on LLPSDB plant entries are not included.
"""

__version__ = "1.4.0"

from .features import FEATURE_NAMES, extract_features
from .model import PlantLLPSModel
from .explain import explain_sequence, logistic_shap, grammar_attributions

__all__ = [
    "FEATURE_NAMES",
    "extract_features",
    "PlantLLPSModel",
    "explain_sequence",
    "logistic_shap",
    "grammar_attributions",
    "__version__",
]
