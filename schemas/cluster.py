from pydantic import BaseModel
from typing import Literal


class ClusterSetting(BaseModel):
    feature: Literal["original", "outcomes", "all"] = "original"
    cluster_method: Literal["rxn_class", "hdbscan", "kmeans"] = "hdbscan"
    fp_type: Literal["morgan"] = "morgan"
    fp_length: int = 512
    fp_radius: int = 1
    classification_threshold: float = 0.2
