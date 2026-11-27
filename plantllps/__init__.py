"""PlantLLPS — plant-aware LLPS propensity scoring (prototype).

Honest scope: this package ships a runnable feature pipeline, a tiny
bundled training set, a transparent composite + logistic stub model,
and a GrapeCondensatome SQLite schema. Full proteome scoring and
production weights trained on LLPSDB plant entries are not included.
"""

__version__ = "0.5.0"

from .features import FEATURE_NAMES, extract_features
from .model import PlantLLPSModel

__all__ = ["FEATURE_NAMES", "extract_features", "PlantLLPSModel", "__version__"]
