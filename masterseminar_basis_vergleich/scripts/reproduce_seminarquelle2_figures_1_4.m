function out = reproduce_seminarquelle2_figures_1_4(mode)
%REPRODUCE_SEMINARQUELLE2_FIGURES_1_4 Recreate selected plane-source plots from Seminarquelle 2.
% Selected model families:
%   - Figure 1: PN with orders 3, 7, 51
%   - Figure 2: HFMn / HFPn with orders 4, 8, 52
%   - Figure 3: PMPn with orders 4, 8, 52
%   - Figure 4: convergence plot for PN, HFMn, HFPn, PMPn
%
% Usage:
%   reproduce_seminarquelle2_figures_1_4()
%   reproduce_seminarquelle2_figures_1_4('smoke')

if nargin < 1 || isempty(mode)
    mode = 'paper';
end

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
addpath(fullfile(rootDir, 'src'));

cfg = mm_default_config();

% Run control.
cfg.io.write_pdf = true;
cfg.io.close_figures = true;
cfg.limiter.type = 'paper';
cfg.timing.enabled = false;
cfg.timing.collect_step_breakdown = false;
cfg.timing.print_summary = false;
cfg.models.families = {'PN', 'HFMn', 'HFPn', 'PMPn'};

% Use a local process pool when the Parallel Computing Toolbox is available.
use_parallel = has_parallel_toolbox();
if use_parallel
    reset_parallel_pool();
    cfg.parallel.enabled = true;
    cfg.parallel.mode = 'single_run';
    cfg.parallel.pool_type = 'local';
    cfg.parallel.num_workers = default_num_workers();
end

[cfg, run_tag, profile_orders] = configure_mode(cfg, mode);
cfg.paths.results = fullfile(rootDir, 'results', ['seminarquelle2_figures_1_4_' run_tag]);
if ~exist(cfg.paths.results, 'dir')
    mkdir(cfg.paths.results);
end

fprintf('[seminarquelle2] mode=%s\n', run_tag);
fprintf('[seminarquelle2] results=%s\n', cfg.paths.results);
fprintf('[seminarquelle2] parallel=%d\n', cfg.parallel.enabled);

res = run_paper2_plane_source(cfg);

fig1_panels = make_panel_spec('PN', profile_orders.PN, '(a) P_N', 'x', '\rho');
fig2_panels = [ ...
    make_panel_spec('HFMn', profile_orders.HFMn, '(a) HFM_n', 'z', 'u_0'), ...
    make_panel_spec('HFPn', profile_orders.HFPn, '(b) HFP_n', 'z', 'u_0')];
fig3_panels = make_panel_spec('PMPn', profile_orders.PMPn, '(a) PMP_n', 'z', 'u_0');

out.figure1 = mm_plot_seminarquelle2_profiles( ...
    cfg.paths.results, fig1_panels, fullfile(cfg.paths.results, 'seminarquelle2_fig1_pn'), cfg.io);
out.figure2 = mm_plot_seminarquelle2_profiles( ...
    cfg.paths.results, fig2_panels, fullfile(cfg.paths.results, 'seminarquelle2_fig2_hf'), cfg.io);
out.figure3 = mm_plot_seminarquelle2_profiles( ...
    cfg.paths.results, fig3_panels, fullfile(cfg.paths.results, 'seminarquelle2_fig3_pmp'), cfg.io);
out.figure4 = mm_plot_seminarquelle2_convergence( ...
    res.table, {'HFMn', 'HFPn', 'PMPn', 'PN'}, ...
    fullfile(cfg.paths.results, 'seminarquelle2_fig4_convergence'), cfg.io);

out.results_dir = cfg.paths.results;
out.csv = res.csv;
out.table = res.table;
out.mode = run_tag;

fprintf('[seminarquelle2] wrote Figure 1-4 to %s\n', cfg.paths.results);

end

function [cfg, run_tag, profile_orders] = configure_mode(cfg, mode)
profile_orders = struct();

switch lower(char(mode))
    case 'paper'
        run_tag = 'paper';
        cfg.paper2.n_cells = 1000;
        cfg.reference.n_cells = 300;
        cfg.reference.n_mu = 64;
        cfg.models.orders.PN = 3:2:51;
        cfg.models.orders.HFMn = 4:2:52;
        cfg.models.orders.HFPn = 4:2:52;
        cfg.models.orders.PMPn = 4:2:52;

        profile_orders.PN = [3, 7, 51];
        profile_orders.HFMn = [4, 8, 52];
        profile_orders.HFPn = [4, 8, 52];
        profile_orders.PMPn = [4, 8, 52];

    case 'smoke'
        run_tag = 'smoke';
        cfg.paper2.n_cells = 120;
        cfg.reference.n_cells = 360;
        cfg.reference.n_mu = 96;
        cfg.models.orders.PN = [3, 7];
        cfg.models.orders.HFMn = [4, 8];
        cfg.models.orders.HFPn = [4, 8];
        cfg.models.orders.PMPn = [4, 8];

        profile_orders.PN = [3, 7];
        profile_orders.HFMn = [4, 8];
        profile_orders.HFPn = [4, 8];
        profile_orders.PMPn = [4, 8];

    otherwise
        error('Unknown mode "%s". Use "paper" or "smoke".', mode);
end

cfg.reference.tf = cfg.paper2.tf;
cfg.reference.domain = cfg.paper2.domain;

end

function panel = make_panel_spec(family, orders, panel_title, x_label, y_label)
panel = struct();
panel.family = char(family);
panel.orders = orders(:).';
panel.panel_title = char(panel_title);
panel.x_label = char(x_label);
panel.y_label = char(y_label);
end

function tf = has_parallel_toolbox()
tf = false;
try
    tf = license('test', 'Distrib_Computing_Toolbox');
catch
end
end

function n = default_num_workers()
n = 8;
try
    n = min(8, feature('numcores'));
catch
end
end

function reset_parallel_pool()
if exist('gcp', 'file') == 2 || exist('gcp', 'builtin') == 5
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
end
end
