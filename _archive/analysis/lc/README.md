# LC analysis archive

`lc_run_pipeline.py` is the earlier in-process runner that mixed LC extraction,
profile building, and plotting in one command. It is retained intact for
reference; the active workflow now uses `preprocessing/preprocessing_run_pipeline.py`,
`analysis/core_analysis_run_pipeline.py`, `analysis/lc/lc_analysis_run_pipeline.py`,
and `plotting/plotting_run_pipeline.py`.
