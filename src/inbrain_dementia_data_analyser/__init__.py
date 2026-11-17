"""inbrain_dementia_data_analyser - harmonised multi-cohort dementia data and study cohorts.

Typical use::

    from inbrain_dementia_data_analyser import (
        enable_logging, validate_all, preprocess_dataset, integrate,
        create_cross_sectional_cohort, create_longitudinal_cohort)

    enable_logging()
    root = "/content/drive/MyDrive/.../Tables_projects_dementia"
    print(validate_all(root)["ADNI"].summary())
    preprocess_dataset("ADNI", root)
    paths = integrate(root)
    cohort = create_longitudinal_cohort(paths["expanded"], modality_label="MRI_T1w")

    from inbrain_dementia_data_analyser.figures import make_cohort_figures
    make_cohort_figures(cohort, "figures", expanded=paths["expanded"])   # needs matplotlib
"""
from ._version import __version__
from .cohorts import (
    CohortResult,
    create_anchor_cohort,
    create_cross_sectional_cohort,
    create_longitudinal_cohort,
    create_longitudinal_multimodal_cohort,
    create_multimodal_cohort,
)
from .config import ConfigError, ProjectConfig, load_config
from .integration import expand_master, integrate, integrated_paths, merge_cohorts
from .logging_utils import enable_logging
from .preprocess import preprocess_all, preprocess_dataset
from .validation import ValidationError, ValidationReport, check_root, validate_all, validate_dataset

__all__ = [
    "__version__", "CohortResult", "ConfigError", "ProjectConfig", "ValidationError",
    "ValidationReport", "check_root", "create_anchor_cohort", "create_cross_sectional_cohort",
    "create_longitudinal_cohort", "create_longitudinal_multimodal_cohort", "create_multimodal_cohort",
    "enable_logging", "expand_master", "integrate", "integrated_paths", "load_config",
    "merge_cohorts", "preprocess_all", "preprocess_dataset", "validate_all", "validate_dataset",
]
