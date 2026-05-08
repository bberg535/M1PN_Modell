function paths = setup_project_paths()
%SETUP_PROJECT_PATHS Configure MATLAB paths for this project layout.

    this_file = mfilename('fullpath');
    if isempty(this_file)
        error('setup_project_paths:noPath', ...
            'Could not resolve the path of setup_project_paths.m.');
    end

    project_root = fileparts(this_file);
    src_dir = fullfile(project_root, 'src');
    scripts_dir = fullfile(project_root, 'scripts');
    results_dir = fullfile(project_root, 'results');

    addpath(src_dir);
    addpath(scripts_dir);

    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    paths = struct();
    paths.project_root = project_root;
    paths.src_dir = src_dir;
    paths.scripts_dir = scripts_dir;
    paths.results_dir = results_dir;
end
