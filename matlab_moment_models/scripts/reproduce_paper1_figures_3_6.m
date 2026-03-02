addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

cfg = mm_default_config();

% Focus profile for Seminarquelle 1 slab tests (Figure 3-6 reconstruction).
cfg.paper1.eval_mu = linspace(-1, 1, 2001);
cfg.io.close_figures = false;
cfg.paper1.figure3_use_cfg_registry = false;

res = run_paper1_figures_3_6(cfg);
disp(res.csv);
