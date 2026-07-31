# archive

`_archive/` holds older scripts outside the current pipeline. The wrappers do not call them; the files remain for reference and may still contain old paths, duplicated helpers, unfinished branches, or imports from the earlier layout.

## folders

- `analysis/`: archived analysis scripts
- `behaviour_code/`: older behaviour scripts
- `exploratory/`: exploratory work and side branches
- `HPC_code/`: older HPC analysis scripts
- `imaging_code/`: older imaging scripts
- `LC_code/`: older LC analysis scripts
- `modelling_code/`: earlier model notebooks
- `preprocessing/`: older preprocessing entry points
- `utils/`: helper modules used by archived scripts

The former in-repository FibreSight workbench is kept under
`preprocessing/imaging/fibre_sight/in_repo_workbench/`. The active workbench and
GUI now live in the separate `dinghaoluo/fibre-sight` repository:
https://github.com/dinghaoluo/fibre-sight

FibreSegger versions 1–5 are also retained under
`preprocessing/imaging/fibre_segger_GUI/` as historical versions of the manual
ROI workflow now handled by the standalone application.
