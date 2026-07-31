function matlab_spike_pipeline_prepare_workdir(filename, raw_recording_dir, output_dir)
    if ~exist(raw_recording_dir, 'dir')
        error('Raw recording folder does not exist: %s', raw_recording_dir);
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    raw_files = dir(raw_recording_dir);
    staged = 0;
    for i = 1:numel(raw_files)
        raw_file = raw_files(i);
        if raw_file.isdir
            continue
        end

        source_file = fullfile(raw_file.folder, raw_file.name);
        [~, ~, ext] = fileparts(raw_file.name);

        % do not seed the output folder with old generated MATLAB products
        if strcmpi(ext, '.mat')
            continue
        end

        % avoid triggering the old dat.orig relinking in the repo output folder
        if endsWith(raw_file.name, '.dat.orig', 'IgnoreCase', true)
            continue
        end

        target_file = fullfile(output_dir, raw_file.name);
        % deleting an old staged link removes only the repo-local entry
        if exist(target_file, 'file')
            delete(target_file);
        end
        % copy rather than link so later processing cannot write into the raw folder
        copyfile(source_file, target_file);
        staged = staged + 1;
    end

    dat_target = fullfile(output_dir, [filename '.dat']);
    dat_orig_source = fullfile(raw_recording_dir, [filename '.dat.orig']);
    if ~exist(dat_target, 'file') && exist(dat_orig_source, 'file')
        copyfile(dat_orig_source, dat_target);
        staged = staged + 1;
    end

    manifest_path = fullfile(output_dir, 'raw_source_path.txt');
    fid = fopen(manifest_path, 'w');
    fprintf(fid, '%s\n', raw_recording_dir);
    fclose(fid);

    disp(['Staged ' num2str(staged) ' raw input files in ' output_dir])
end
