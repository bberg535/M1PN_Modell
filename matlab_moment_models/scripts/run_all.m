function run_all()
%RUN_ALL Execute both benchmark scripts with default configuration.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'src'));

cfg = mm_default_config();
%cfg.limiter.type = 'mcl';
cfg.parallel.enabled = true;
cfg.parallel.mode = 'single-run';
cfg.parallel.pool_type = 'local';
cfg.parallel.num_workers = 8;


%disp('Running paper 1 (Gauss/Heaviside)...');
%r1 = run_paper1_gauss_heaviside(cfg);
%disp(r1.table);

disp('Running paper 2 (Plane source)...');
r2 = run_paper2_plane_source(cfg);
disp(r2.table);

end
