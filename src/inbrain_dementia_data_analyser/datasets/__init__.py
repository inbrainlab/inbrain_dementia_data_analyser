"""One pipeline module per dataset; each exposes ``run(ctx) -> DataFrame``."""
from .base import DatasetContext, dataset_module, run_dataset

__all__ = ["DatasetContext", "dataset_module", "run_dataset"]
