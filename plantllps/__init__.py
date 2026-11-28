"""PlantLLPS — plant-aware LLPS propensity scoring (prototype).

Honest scope: this package ships a runnable feature pipeline, a tiny
bundled training set, a transparent composite + logistic stub model,
and a GrapeCondensatome SQLite schema. Full proteome scoring and
production weights trained on LLPSDB plant entries are not included.
"""

from pathlib import Path

_VERSION_FILE = Path(__file__).resolve().parent.parent / "VERSION"
try:
    __version__ = _VERSION_FILE.read_text(encoding="utf-8").strip() or "1.11.0"
except OSError:
    __version__ = "1.11.0"

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
