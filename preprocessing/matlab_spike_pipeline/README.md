# MATLAB spike pipeline

`runspikepipeline.m` and `runspikepipeline_pix.m` are the two entry points I added on 23 January 2023. They call the working Wang lab `support/` tree, where older code by Yingxue Wang, Dinghao Luo, and other lab contributors sits alongside named utilities such as the Circular Statistics Toolbox for MATLAB. The support tree still carries its original mixture of spelling, local paths, backup folders, and dated edits.

The two entry points now resolve the recording name, staging directory and generated output paths directly, then write generated output under `data/matlab_spike_pipeline/`. The five short `support/matlab_spike_pipeline_*.m` helpers that only resolved a name or path, or copied one file, have gone; the equally thin `pathmatlab_smTr_dinghao.m` wrapper has gone too, so these steps can be read where the pipeline is actually worked through.
