function cmp = compare_limiters_mcl_vs_paper(cfg)
%COMPARE_LIMITERS_MCL_VS_PAPER Compare paper limiter against MCL limiter.

if nargin < 1 || isempty(cfg)
    cfg = mm_default_config();
end

if ~exist(cfg.paths.results, 'dir')
    mkdir(cfg.paths.results);
end

closeFigures = get_close_figures(cfg);

cfgPaper = cfg;
cfgPaper.limiter.type = 'paper';
cfgPaper.limiter.paper_lp_characteristic = true;

cfgMCL = cfg;
cfgMCL.limiter.type = 'mcl';

resPaper = run_paper2_plane_source(cfgPaper);
resMCL = run_paper2_plane_source(cfgMCL);

Tp = resPaper.table;
Tm = resMCL.table;

Tp = renamevars(Tp, {'L1', 'Linf', 'symmetry_L1_rel', 'min_rho', 'runtime_s', 'n_steps'}, ...
    {'L1_paper', 'Linf_paper', 'sym_paper', 'minrho_paper', 'runtime_paper', 'steps_paper'});
Tm = renamevars(Tm, {'L1', 'Linf', 'symmetry_L1_rel', 'min_rho', 'runtime_s', 'n_steps'}, ...
    {'L1_mcl', 'Linf_mcl', 'sym_mcl', 'minrho_mcl', 'runtime_mcl', 'steps_mcl'});

Tcmp = innerjoin(Tp, Tm, 'Keys', {'model', 'order', 'nMom'});
Tcmp.dL1_abs = Tcmp.L1_mcl - Tcmp.L1_paper;
Tcmp.dL1_rel = Tcmp.dL1_abs ./ max(Tcmp.L1_paper, eps);
Tcmp.dLinf_abs = Tcmp.Linf_mcl - Tcmp.Linf_paper;
Tcmp.dRuntime_abs = Tcmp.runtime_mcl - Tcmp.runtime_paper;

outCsv = fullfile(cfg.paths.results, 'limiter_comparison_paper_vs_mcl.csv');
writetable(Tcmp, outCsv);

plot_comparison(Tcmp, cfg.paths.results, closeFigures);

cmp = struct();
cmp.table = Tcmp;
cmp.csv = outCsv;
cmp.paper = resPaper;
cmp.mcl = resMCL;

end

function plot_comparison(Tcmp, outDir, closeFigures)
fig = figure('Color', 'w');
tiledlayout(1, 2);

labels = strcat(string(Tcmp.model), "-", string(Tcmp.order));

nexttile;
b = bar(categorical(labels), [Tcmp.L1_paper, Tcmp.L1_mcl], 'grouped');
b(1).DisplayName = 'paper';
b(2).DisplayName = 'mcl';
ylabel('L1 error');
title('Plane source: L1 comparison');
grid(gca, 'on');
legend('Location', 'best');

nexttile;
b2 = bar(categorical(labels), [Tcmp.runtime_paper, Tcmp.runtime_mcl], 'grouped');
b2(1).DisplayName = 'paper';
b2(2).DisplayName = 'mcl';
ylabel('runtime [s]');
title('Plane source: runtime comparison');
grid(gca, 'on');
legend('Location', 'best');

outBase = fullfile(outDir, 'limiter_comparison_paper_vs_mcl');
saveas(fig, [outBase '.png']);
try
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
catch
end
if closeFigures
    close(fig);
end
end

function closeFigures = get_close_figures(cfg)
closeFigures = true;
if isfield(cfg, 'io') && isfield(cfg.io, 'close_figures')
    closeFigures = logical(cfg.io.close_figures);
end
end
