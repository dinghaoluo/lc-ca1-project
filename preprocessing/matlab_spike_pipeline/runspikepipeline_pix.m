function runspikepipeline_pix()
    % Runs all pixel-rig spike-sorting pipeline functions.
    % Dinghao 23 Jan, 2023

    pipeline_dir = fileparts(mfilename('fullpath'));
    support_dir = fullfile(pipeline_dir, 'support');
    addpath(genpath(support_dir));
    addpath(fullfile(support_dir, 'recording_lists'), '-begin');
    original_dir = pwd;
    cleanup_obj = onCleanup(@() cd(original_dir));
    output_parent_dir = fullfile(fileparts(fileparts(pipeline_dir)), ...
        'data', 'matlab_spike_pipeline');

    listRecordingsActiveLickPathHPCLCOpt = [];
    listRecordingsActiveLickFileNameHPCLCOpt = [];
    RecordingList_dinghao;
    paths = listRecordingsActiveLickPathHPCLCOpt(1:end, :);
    files = listRecordingsActiveLickFileNameHPCLCOpt(1:end, :);
    tot_rec = size(paths, 1);

    for i = 1:tot_rec
        pathname = strtrim(paths(i, :));
        filename = strtrim(files(i, :));
        suffix_idx = strfind(filename, '_DataStructure');
        if ~isempty(suffix_idx)
            filename = filename(1:suffix_idx(1)-1);
        end

        raw_recording_dir = strtrim(char(pathname));
        while ~isempty(raw_recording_dir) && ...
                (raw_recording_dir(end) == filesep || ...
                 raw_recording_dir(end) == '/' || raw_recording_dir(end) == '\')
            raw_recording_dir = raw_recording_dir(1:end-1);
        end
        output_dir = fullfile(output_parent_dir, filename);

        disp(filename);
        disp(['Raw recording path ' raw_recording_dir])
        disp(['Writing outputs to ' output_dir])

        matlab_spike_pipeline_prepare_workdir(filename, raw_recording_dir, output_dir);
        cd(output_dir);

        fullfilename = [filename '_DataStructure_mazeSection1_TrialType1'];

        % generate B.mat, BTDT.mat, -whl.mat, .NeuronQuality.mat,
        % _BehavElectrDataLFP.mat, _BehavElectrDataLFP_CCG.mat, _eeg_1250Hz.mat.
        GenerateBehavElectroDataStructures_smTrMPFIv2_opto_pix(output_parent_dir, filename);

        % Generate _DataStructure_mazeSection1_TrialType1.mat.
        % DOUBLE CHECK STIM PULSE NUMBER.
        GetTrials_smTrMPFI(filename, 1, 1, 0, 0);

        % generate processed Info, ext, im, runSpeed, FR_Run0, Depth,
        % convSpikesTime9p6ms, convSpikesDist20mm, PeakFR20mm, ThetaPhaseH,
        % ThetaPhaseL, Concatsp, CCG, ThetaMod, SpInfo, FieldWidthLR_20mm_L,
        % burstAll_THH, burstAll_THL.
        disp([newline 'PROCESSING'])
        ProcessingMice_smTr('./', fullfilename, 0);

        % generate processed thetaPower, alignRun_msess, alignRew_msess,
        % alignCue_msess, alignCueOff_msess, alignedSpikesPerNPerT_msess,
        % convSpikesAligned_msess, PeakFRAligned_msess,
        % thetaPhaseOverTimeligned_msess, lickDist_msess, runSpeedDist_msess,
        % FRAlignedRun_msess, ThetaPhaseAligned files, and rasters.
        disp([newline 'ALIGNING'])
        ProcessingAlignedWithStim('./', fullfilename, 0, 1, 0, 3);

        disp([newline 'PROCESSING CONTROL ONLY'])
        ProcessingMice_smTrCtrlOnly('./', fullfilename, 0, 1)

        % tagging-related functions are disabled for pixel recordings
    %     stimEffect_NewData_MPFI('./', fullfilename);
    %     plotStimRasterWrapper('./', fullfilename, -1);
    end
end
