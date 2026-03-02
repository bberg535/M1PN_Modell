function run_all()
%RUN_ALL Execute both benchmark scripts with default configuration.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'src'));

cfg = mm_default_config();

disp('Running paper 1 (Gauss/Heaviside)...');
r1 = run_paper1_gauss_heaviside(cfg);
disp(r1.table);

disp('Running paper 2 (Plane source)...');
r2 = run_paper2_plane_source(cfg);
disp(r2.table);

end
