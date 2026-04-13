% REPRODUCE_SEMINARQUELLE2_FIGURES_1_4
% Recreates Seminarquelle 2, Figure 1-4 (plane-source case).
% Toggle very high orders via include_very_high_orders.

include_very_high_orders = false;   % true: include M/P 51 and H*/PM* 52 and extended convergence orders
keep_figures_open = true;           % false: close figures after saving
write_pdf = true;                   % also export vector PDF
limiter_type = 'paper';             % 'paper' or 'mcl'

this_dir = fileparts(mfilename('fullpath'));
mm_root = fullfile(this_dir, 'matlab_moment_models');
addpath(fullfile(mm_root, 'src'));

cfg0 = mm_default_config();
cfg0.io.close_figures = ~logical(keep_figures_open);
cfg0.io.write_pdf = logical(write_pdf);
cfg0.limiter.type = limiter_type;
cfg0.timing.enabled = false;
cfg0.timing.collect_step_breakdown = false;
cfg0.timing.print_summary = false;
cfg0.timing.print_progress = false;

cfg0.paper2.domain = [-1.2, 1.2];
cfg0.paper2.tf = 1.0;
cfg0.paper2.n_cells = 80;
cfg0.reference.n_cells = 200;
cfg0.reference.n_mu = 64;
cfg0.reference.tf = cfg0.paper2.tf;
cfg0.reference.domain = cfg0.paper2.domain;

orders = build_orders(include_very_high_orders);
mode_tag = mode_string(include_very_high_orders);
out_dir = fullfile(mm_root, 'results', ['seminarquelle2_fig1_4_' mode_tag]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('Output directory: %s\n', out_dir);

cfg_fig = cfg0;
cfg_fig.paths.results = out_dir;

% Figure 1: MN / PN profiles.
cfg1 = cfg_fig;
cfg1.models.families = {'MN', 'PN'};
cfg1.models.orders.MN = orders.fig1_mn;
cfg1.models.orders.PN = orders.fig1_pn;
run_paper2_plane_source(cfg1);
plot_pair_figure(cfg1.paths.results, 'MN', orders.fig1_mn, 'PN', orders.fig1_pn, ...
    'Seminarquelle 2 - Figure 1', fullfile(cfg1.paths.results, ['seminarquelle2_fig1_' mode_tag]), cfg1.io);

% Figure 2: HFMn / HFPn profiles.
cfg2 = cfg_fig;
cfg2.models.families = {'HFMn', 'HFPn'};
cfg2.models.orders.HFMn = orders.fig2_hfm;
cfg2.models.orders.HFPn = orders.fig2_hfp;
run_paper2_plane_source(cfg2);
plot_pair_figure(cfg2.paths.results, 'HFMn', orders.fig2_hfm, 'HFPn', orders.fig2_hfp, ...
    'Seminarquelle 2 - Figure 2', fullfile(cfg2.paths.results, ['seminarquelle2_fig2_' mode_tag]), cfg2.io);

% Figure 3: PMMn / PMPn profiles.
cfg3 = cfg_fig;
cfg3.models.families = {'PMMn', 'PMPn'};
cfg3.models.orders.PMMn = orders.fig3_pmm;
cfg3.models.orders.PMPn = orders.fig3_pmp;
run_paper2_plane_source(cfg3);
plot_pair_figure(cfg3.paths.results, 'PMMn', orders.fig3_pmm, 'PMPn', orders.fig3_pmp, ...
    'Seminarquelle 2 - Figure 3', fullfile(cfg3.paths.results, ['seminarquelle2_fig3_' mode_tag]), cfg3.io);

% Figure 4: convergence L1 / Linf.
cfg4 = cfg_fig;
cfg4.models.families = {'HFMn', 'HFPn', 'PMMn', 'PMPn', 'MN', 'PN'};
cfg4.models.orders.HFMn = orders.conv_hfm;
cfg4.models.orders.HFPn = orders.conv_hfp;
cfg4.models.orders.PMMn = orders.conv_pmm;
cfg4.models.orders.PMPn = orders.conv_pmp;
cfg4.models.orders.MN = orders.conv_mn;
cfg4.models.orders.PN = orders.conv_pn;
res4 = run_paper2_plane_source(cfg4);
plot_convergence_figure(res4.table, ...
    {'HFMn', 'HFPn', 'PMMn', 'PMPn', 'MN', 'PN'}, ...
    'Seminarquelle 2 - Figure 4', fullfile(cfg4.paths.results, ['seminarquelle2_fig4_' mode_tag]), cfg4.io);

fprintf('Finished: Figure 1-4 written to %s\n', out_dir);

function orders = build_orders(include_high)
orders = struct();
if include_high
    orders.fig1_mn = [3, 7, 51];
    orders.fig1_pn = [3, 7, 51];
    orders.fig2_hfm = [4, 8, 52];
    orders.fig2_hfp = [4, 8, 52];
    orders.fig3_pmm = [4, 8, 52];
    orders.fig3_pmp = [4, 8, 52];

    orders.conv_mn = 3:2:51;
    orders.conv_pn = 3:2:51;
    orders.conv_hfm = 4:2:52;
    orders.conv_hfp = 4:2:52;
    orders.conv_pmm = 4:2:52;
    orders.conv_pmp = 4:2:52;
else
    orders.fig1_mn = [3, 7];
    orders.fig1_pn = [3, 7];
    orders.fig2_hfm = [4, 8];
    orders.fig2_hfp = [4, 8];
    orders.fig3_pmm = [4, 8];
    orders.fig3_pmp = [4, 8];

    orders.conv_mn = [3, 7];
    orders.conv_pn = [3, 7];
    orders.conv_hfm = [4, 8];
    orders.conv_hfp = [4, 8];
    orders.conv_pmm = [4, 8];
    orders.conv_pmp = [4, 8];
end
end

function s = mode_string(include_high)
if include_high
    s = 'high';
else
    s = 'base';
end
end

function plot_pair_figure(results_dir, left_model, left_orders, right_model, right_orders, fig_title, out_base, io_cfg)
fig = figure('Color', 'w', 'Name', fig_title);
tiledlayout(1, 2);

nexttile;
plot_profiles(results_dir, left_model, left_orders);
title(['(a) ' family_display_name(left_model)]);

nexttile;
plot_profiles(results_dir, right_model, right_orders);
title(['(b) ' family_display_name(right_model)]);

saveas(fig, [out_base '.png']);
if io_cfg.write_pdf
    try
        exportgraphics(fig, [out_base '.pdf'], 'ContentType', 'vector');
    catch
    end
end
if io_cfg.close_figures
    close(fig);
end
end

function plot_profiles(results_dir, model_name, orders)
hold on;
grid(gca, 'on');
ref_plotted = false;

for k = 1:numel(orders)
    ord = orders(k);
    prof = read_profile(results_dir, model_name, ord);
    mask = prof.z >= 0;

    if ~ref_plotted
        plot(prof.z(mask), prof.rho_ref(mask), 'k--', 'LineWidth', 1.4, 'DisplayName', 'Ref');
        ref_plotted = true;
    end

    plot(prof.z(mask), prof.rho_model(mask), 'LineWidth', 1.2, ...
        'DisplayName', model_curve_label(model_name, ord));
end

xlabel('z');
ylabel('rho');
xlim([0, 1]);
legend('Location', 'best');
end

function T = read_profile(results_dir, model_name, order)
fn = fullfile(results_dir, sprintf('paper2_plane_source_%s_%d_profile.csv', model_name, order));
if ~exist(fn, 'file')
    error('Missing profile file: %s', fn);
end
T = readtable(fn);
end

function plot_convergence_figure(T, model_order, fig_title, out_base, io_cfg)
fig = figure('Color', 'w', 'Name', fig_title);
tiledlayout(1, 2);

nexttile;
hold on;
for i = 1:numel(model_order)
    m = model_order{i};
    idx = strcmp(string(T.model), string(m));
    Ti = T(idx, :);
    if isempty(Ti)
        continue;
    end
    [x, p] = sort(Ti.nMom);
    y = Ti.L1(p);
    loglog(x, y, '-o', 'LineWidth', 1.2, 'DisplayName', family_display_name(m));
end
grid(gca, 'on');
xlabel('n');
ylabel('L1 error');
title('(a) L1 convergence');
legend('Location', 'best');

nexttile;
hold on;
for i = 1:numel(model_order)
    m = model_order{i};
    idx = strcmp(string(T.model), string(m));
    Ti = T(idx, :);
    if isempty(Ti)
        continue;
    end
    [x, p] = sort(Ti.nMom);
    y = Ti.Linf(p);
    loglog(x, y, '-o', 'LineWidth', 1.2, 'DisplayName', family_display_name(m));
end
grid(gca, 'on');
xlabel('n');
ylabel('Linf error');
title('(b) Linf convergence');
legend('Location', 'best');

saveas(fig, [out_base '.png']);
if io_cfg.write_pdf
    try
        exportgraphics(fig, [out_base '.pdf'], 'ContentType', 'vector');
    catch
    end
end
if io_cfg.close_figures
    close(fig);
end
end

function s = family_display_name(model_name)
s = char(model_name);
if strcmp(s, 'MN')
    s = 'M_N';
elseif strcmp(s, 'PN')
    s = 'P_N';
end
end

function s = model_curve_label(model_name, order)
base = char(model_name);
if strcmp(base, 'MN')
    s = sprintf('M%d', order);
elseif strcmp(base, 'PN')
    s = sprintf('P%d', order);
elseif strcmp(base, 'HFMn')
    s = sprintf('HFM%d', order);
elseif strcmp(base, 'HFPn')
    s = sprintf('HFP%d', order);
elseif strcmp(base, 'PMMn')
    s = sprintf('PMM%d', order);
elseif strcmp(base, 'PMPn')
    s = sprintf('PMP%d', order);
else
    s = sprintf('%s%d', base, order);
end
end
