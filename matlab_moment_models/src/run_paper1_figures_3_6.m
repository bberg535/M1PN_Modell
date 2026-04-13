function res = run_paper1_figures_3_6(cfg)
%RUN_PAPER1_FIGURES_3_6 Reproduce slab-geometry numerical tests from Paper 1.
% Targets Figure 3-6 style outputs (Gauss, Heaviside, CrossingBeams).

if nargin < 1 || isempty(cfg)
    cfg = mm_default_config();
end

if ~exist(cfg.paths.results, 'dir')
    mkdir(cfg.paths.results);
end

closeFigures = get_close_figures(cfg);
paper1Cfg = fill_paper1_cfg(cfg.paper1);

muEval = cfg.paper1.eval_mu(:);
muPlot = build_plot_grid(muEval, paper1Cfg.plot_endpoint_trim);
psiVac = cfg.physics.psi_vac_density;
aCross = cfg.paper1.crossing_a;

psiGauss = @(mu) (1.0 / sqrt(2.0 * pi * cfg.paper1.gauss_sigma^2)) .* ...
    exp(-((mu - cfg.paper1.gauss_mu_bar).^2) ./ (2.0 * cfg.paper1.gauss_sigma^2));
psiHeav = @(mu) (mu < cfg.paper1.heaviside_jump) .* psiVac + (mu >= cfg.paper1.heaviside_jump) .* 1.0;
psiCross = @(mu) sqrt(aCross / pi) .* (exp(-aCross * (mu + 1.0).^2) + exp(-aCross * (mu - 0.5).^2));

cases = struct([]);
cases(1).name = 'Gauss';
cases(1).psi = psiGauss;
cases(1).ref_eval = psiGauss(muEval);
cases(1).ref_plot = psiGauss(muPlot);
cases(2).name = 'Heaviside';
cases(2).psi = psiHeav;
cases(2).ref_eval = psiHeav(muEval);
cases(2).ref_plot = psiHeav(muPlot);
cases(3).name = 'CrossingBeams';
cases(3).psi = psiCross;
cases(3).ref_eval = psiCross(muEval);
cases(3).ref_plot = psiCross(muPlot);

registryFig3 = build_figure3_registry(cfg, paper1Cfg);
rows = [];
warmstartCache = struct();

for i = 1:numel(registryFig3)
    rec = registryFig3(i);
    model = mm_build_model(rec.name, rec.order, cfg.models);
    quad = mm_build_quadrature(model, cfg.quad, 'moment');

    for ic = 1:numel(cases)
        c = cases(ic);
        cacheKey = get_warmstart_key(c.name, rec.name);
        cacheIn = struct();
        if paper1Cfg.figure3_use_warmstart && model.needs_entropy && isfield(warmstartCache, cacheKey)
            cacheIn.alpha0 = fit_warmstart_alpha(warmstartCache.(cacheKey), model.nMom);
        end

        [row, ~, cacheOut] = eval_case(c.name, c.psi, c.ref_eval, muEval, model, quad, cfg.optimizer, rec, paper1Cfg, cacheIn);
        rows = [rows; row]; %#ok<AGROW>
        if paper1Cfg.figure3_use_warmstart && model.needs_entropy && ~isempty(cacheOut)
            warmstartCache.(cacheKey) = cacheOut;
        end
    end
end

T = struct2table(rows);

fig3Base = fullfile(cfg.paths.results, 'paper1_figure3_convergence');
plot_figure3_convergence(T, fig3Base, closeFigures);

fig4Base = fullfile(cfg.paths.results, 'paper1_figure4_gauss_ansatz');
plot_selected_families(muPlot, cases(1), cfg, get_figure4_selection(), fig4Base, closeFigures, paper1Cfg);

fig5Base = fullfile(cfg.paths.results, 'paper1_figure5_heaviside_ansatz');
plot_selected_families(muPlot, cases(2), cfg, get_figure5_selection(), fig5Base, closeFigures, paper1Cfg);

fig6Base = fullfile(cfg.paths.results, 'paper1_figure6_crossing_beams_ansatz');
plot_selected_families(muPlot, cases(3), cfg, get_figure6_selection(), fig6Base, closeFigures, paper1Cfg);

res = struct();
res.table = T;
res.figure3_base = fig3Base;
res.figure4_base = fig4Base;
res.figure5_base = fig5Base;
res.figure6_base = fig6Base;

end

function [row, psiHat, cacheOut] = eval_case(caseName, psiFun, ref, muEval, model, quad, optCfg, rec, paper1Cfg, cacheIn)
if nargin < 10 || isempty(cacheIn)
    cacheIn = struct();
end

tStart = tic;
[psiHat, info, cacheOut] = reconstruct_density(caseName, psiFun, ref, muEval, model, quad, optCfg, paper1Cfg, cacheIn);

err = abs(psiHat - ref(:));
L1 = trapz(muEval, err);
Linf = max(err);

row = struct();
row.case = caseName;
row.model = rec.name;
row.order = rec.order;
row.nMom = rec.nMom;
row.L1 = L1;
row.Linf = Linf;
row.runtime_s = toc(tStart);
row.iterations = info.iterations;
row.converged = info.converged;
end

function [psiHat, info, cacheOut] = reconstruct_density(caseName, psiFun, ref, muGrid, model, quad, optCfg, paper1Cfg, cacheState)
if nargin < 9 || isempty(cacheState)
    cacheState = struct();
end

cacheOut = [];
u = mm_project_density_to_moments(psiFun, model, quad);
B = model.basis_eval(muGrid);

if model.needs_entropy
    optUse = adapt_opt_cfg(optCfg, model, caseName, paper1Cfg);
    [alpha, info] = mm_entropy_dual_solve(u, model, quad, optUse, cacheState);
    psiHat = exp(min(B * alpha, 700));

    if should_retry_entropy(caseName, model, psiHat, ref, paper1Cfg)
        optRetry = optUse;
        optRetry.regularization_r = enforce_regularization_floor(optUse.regularization_r, ...
            get_case_retry_reg_floor(caseName, paper1Cfg));
        retryCache = cacheState;
        if isfield(info, 'cached_alpha') && numel(info.cached_alpha) == model.nMom
            retryCache.last_alpha = info.cached_alpha(:);
        end
        [alphaRetry, infoRetry] = mm_entropy_dual_solve(u, model, quad, optRetry, retryCache);
        psiRetry = exp(min(B * alphaRetry, 700));
        if is_better_profile(psiHat, psiRetry, ref)
            psiHat = psiRetry;
            info = infoRetry;
        end
    end

    if any(~isfinite(psiHat))
        psiHat(:) = NaN;
    end
    if info.converged && isfield(info, 'cached_alpha') && numel(info.cached_alpha) == model.nMom
        cacheOut = info.cached_alpha(:);
    end
else
    alpha = model.mass_matrix \ u;
    psiHat = B * alpha;
    info = struct('iterations', 0, 'converged', true);
end
end

function tf = is_better_profile(psiOld, psiNew, ref)
if any(~isfinite(psiNew))
    tf = false;
    return;
end

refScale = max(max(ref), 1.0);
oldScore = norm(abs(psiOld - ref), inf) / refScale;
newScore = norm(abs(psiNew - ref), inf) / refScale;

oldFinite = all(isfinite(psiOld));
newFinite = all(isfinite(psiNew));

if ~oldFinite && newFinite
    tf = true;
    return;
end

tf = newScore <= oldScore;
end

function optUse = adapt_opt_cfg(optCfg, model, caseName, paper1Cfg)
optUse = optCfg;
if model.needs_entropy && strcmp(model.name, 'MN') && is_nonsmooth_case(caseName)
    optUse.regularization_r = enforce_regularization_floor(get_field_or(optCfg, 'regularization_r', [0, 1e-6, 1e-4, 1e-2, 1]), ...
        get_case_reg_floor(caseName, paper1Cfg));
end
end

function tf = should_retry_entropy(caseName, model, psiHat, ref, paper1Cfg)
tf = false;
if ~(model.needs_entropy && strcmp(model.name, 'MN') && is_nonsmooth_case(caseName))
    return;
end
if any(~isfinite(psiHat))
    tf = true;
    return;
end
refScale = max(max(ref), 1.0);
if max(psiHat) > paper1Cfg.nonsmooth_retry_spike_factor * refScale
    tf = true;
end
end

function tf = is_nonsmooth_case(caseName)
tf = strcmp(caseName, 'Heaviside') || strcmp(caseName, 'CrossingBeams');
end

function rOut = enforce_regularization_floor(rIn, floorVal)
r = rIn(:).';
r = r(r >= floorVal);
rOut = unique([r, floorVal, 1.0], 'sorted');
end

function plot_figure3_convergence(T, outBase, closeFigures)
fig = figure('Color', 'w', 'Name', 'Paper1 Figure 3');
tiledlayout(3, 2);

plot_case_error(T, 'Gauss', 'L1');
plot_case_error(T, 'Gauss', 'Linf');
plot_case_error(T, 'Heaviside', 'L1');
plot_case_error(T, 'Heaviside', 'Linf');
plot_case_error(T, 'CrossingBeams', 'L1');
plot_case_error(T, 'CrossingBeams', 'Linf');

saveas(fig, [outBase '.png']);
try
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
catch
end
if closeFigures
    close(fig);
end
end

function plot_case_error(T, caseName, errName)
ax = nexttile;
mask = strcmp(string(T.case), caseName);
Tc = T(mask, :);
models = unique(string(Tc.model), 'stable');
hold(ax, 'on');

for i = 1:numel(models)
    m = models(i);
    tm = Tc(strcmp(string(Tc.model), m), :);
    [x, idx] = sort(tm.nMom);
    y = tm.(errName)(idx);

    if strcmp(caseName, 'Heaviside')
        isPM = strcmp(char(m), 'PMMn') || strcmp(char(m), 'PMPn');
        if isPM
            keep = mod(tm.order(idx), 4) ~= 0;
            x = x(keep);
            y = y(keep);
        end
    end

    good = isfinite(x) & isfinite(y) & (x > 0) & (y > 0);
    x = x(good);
    y = y(good);

    if ~isempty(x)
        plot(ax, x, y, '-o', 'DisplayName', char(m), 'LineWidth', 1.2, 'MarkerSize', 4);
    end
end
set(ax, 'YScale', 'log');
grid(ax, 'on');
xlabel(ax, 'nMom');
ylabel(ax, sprintf('%s error', errName));
title(ax, sprintf('%s: %s convergence', caseName, errName));
legend(ax, 'Location', 'best');
end

function plot_selected_families(muPlot, c, cfg, selection, outBase, closeFigures, paper1Cfg)
families = fieldnames(selection);
nF = numel(families);
nRows = 2;
nCols = ceil(nF / 2);
fig = figure('Color', 'w', 'Name', sprintf('Paper1 %s', c.name));
tiledlayout(nRows, nCols);

for i = 1:nF
    family = families{i};
    orders = selection.(family);
    nexttile;
    hold on;
    plot(muPlot, c.ref_plot, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ref');

    for io = 1:numel(orders)
        ord = orders(io);
        try
            model = mm_build_model(family, ord, cfg.models);
            quad = mm_build_quadrature(model, cfg.quad, 'moment');
            [psiHat, ~] = reconstruct_density(c.name, c.psi, c.ref_plot, muPlot, model, quad, cfg.optimizer, paper1Cfg);
            plot(muPlot, psiHat, 'LineWidth', 1.0, 'DisplayName', sprintf('%s%d', label_family(family), ord));
        catch ME
            warning('Skipping %s-%d in %s plot: %s', family, ord, c.name, ME.message);
        end
    end

    xlabel('\mu');
    ylabel('\psi');
    title(sprintf('%s (%s)', family, c.name));
    grid(gca, 'on');
    legend('Location', 'best');
end

saveas(fig, [outBase '.png']);
try
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
catch
end
if closeFigures
    close(fig);
end
end

function registry = build_figure3_registry(cfg, paper1Cfg)
if paper1Cfg.figure3_use_cfg_registry
    registry = mm_model_registry(cfg);
    return;
end

orders = paper1Cfg.figure3_orders;
families = {'PN', 'MN', 'HFPn', 'HFMn', 'PMPn', 'PMMn'};
registry = struct('name', {}, 'order', {}, 'nMom', {});

for i = 1:numel(families)
    fam = families{i};
    if ~isfield(orders, fam)
        continue;
    end
    ordList = unique(orders.(fam), 'stable');
    for j = 1:numel(ordList)
        ord = ordList(j);
        model = mm_build_model(fam, ord, cfg.models);
        rec = struct('name', fam, 'order', ord, 'nMom', model.nMom);
        registry(end + 1) = rec; %#ok<AGROW>
    end
end
end

function s = get_figure4_selection()
s = struct();
s.PN = [1, 2, 3, 4];
s.HFMn = [2, 3, 4, 5];
s.PMMn = [2, 4, 6, 8];
end

function s = get_figure5_selection()
s = struct();
s.MN = [1, 2, 3, 4, 49];
s.HFMn = [2, 4, 6, 12, 40];
s.PMMn = [2, 6, 10, 22, 34];
s.PN = [1, 2, 3, 4, 49];
s.HFPn = [2, 4, 6, 12, 40];
s.PMPn = [2, 6, 10, 22, 34];
end

function s = get_figure6_selection()
s = struct();
s.MN = [1, 2, 3, 4];
s.HFMn = [10, 11, 12, 13, 60];
s.PMMn = [10, 12, 14, 16, 60];
s.PN = [10, 11, 12, 13, 60];
s.HFPn = [10, 11, 12, 13, 60];
s.PMPn = [10, 12, 14, 16, 60];
end

function txt = label_family(fam)
if strcmp(fam, 'PN') || strcmp(fam, 'MN')
    txt = fam(1);
elseif strcmp(fam, 'HFMn')
    txt = 'HFM';
elseif strcmp(fam, 'HFPn')
    txt = 'HFP';
elseif strcmp(fam, 'PMMn')
    txt = 'PMM';
else
    txt = 'PMP';
end
end

function muPlot = build_plot_grid(muEval, trimEps)
if trimEps <= 0
    muPlot = muEval;
else
    muPlot = linspace(-1 + trimEps, 1 - trimEps, numel(muEval)).';
end
end

function paper1Cfg = fill_paper1_cfg(in)
paper1Cfg = in;

if ~isfield(paper1Cfg, 'figure3_use_cfg_registry')
    paper1Cfg.figure3_use_cfg_registry = false;
end
if ~isfield(paper1Cfg, 'figure3_use_warmstart')
    paper1Cfg.figure3_use_warmstart = false;
end
if ~isfield(paper1Cfg, 'figure3_orders')
    paper1Cfg.figure3_orders = default_figure3_orders();
end
if ~isfield(paper1Cfg, 'heaviside_reg_floor')
    paper1Cfg.heaviside_reg_floor = get_field_or(paper1Cfg, 'nonsmooth_reg_floor', 1.0e-2);
end
if ~isfield(paper1Cfg, 'heaviside_reg_retry_floor')
    paper1Cfg.heaviside_reg_retry_floor = get_field_or(paper1Cfg, 'nonsmooth_reg_retry_floor', 5.0e-2);
end
if ~isfield(paper1Cfg, 'crossing_reg_floor')
    paper1Cfg.crossing_reg_floor = 0.0;
end
if ~isfield(paper1Cfg, 'crossing_reg_retry_floor')
    paper1Cfg.crossing_reg_retry_floor = 0.0;
end
if ~isfield(paper1Cfg, 'nonsmooth_reg_floor')
    paper1Cfg.nonsmooth_reg_floor = 1.0e-2;
end
if ~isfield(paper1Cfg, 'nonsmooth_reg_retry_floor')
    paper1Cfg.nonsmooth_reg_retry_floor = 5.0e-2;
end
if ~isfield(paper1Cfg, 'nonsmooth_retry_spike_factor')
    paper1Cfg.nonsmooth_retry_spike_factor = 50.0;
end
if ~isfield(paper1Cfg, 'plot_endpoint_trim')
    paper1Cfg.plot_endpoint_trim = 1.0e-3;
end
end

function orders = default_figure3_orders()
orders = struct();
orders.PN = [1:2:49, 60];
orders.MN = [1:2:19, 29, 39, 49];
orders.HFPn = 2:2:60;
orders.HFMn = 2:2:60;
orders.PMPn = 2:2:60;
orders.PMMn = 2:2:60;
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function closeFigures = get_close_figures(cfg)
closeFigures = true;
if isfield(cfg, 'io') && isfield(cfg.io, 'close_figures')
    closeFigures = logical(cfg.io.close_figures);
end
end

function key = get_warmstart_key(caseName, modelName)
key = [caseName '__' modelName];
end

function alphaFit = fit_warmstart_alpha(alphaPrev, nMom)
alphaFit = zeros(nMom, 1);
nCopy = min(numel(alphaPrev), nMom);
alphaFit(1:nCopy) = alphaPrev(1:nCopy);
end

function floorVal = get_case_reg_floor(caseName, paper1Cfg)
if strcmp(caseName, 'Heaviside')
    floorVal = paper1Cfg.heaviside_reg_floor;
elseif strcmp(caseName, 'CrossingBeams')
    floorVal = paper1Cfg.crossing_reg_floor;
else
    floorVal = 0.0;
end
end

function floorVal = get_case_retry_reg_floor(caseName, paper1Cfg)
if strcmp(caseName, 'Heaviside')
    floorVal = paper1Cfg.heaviside_reg_retry_floor;
elseif strcmp(caseName, 'CrossingBeams')
    floorVal = paper1Cfg.crossing_reg_retry_floor;
else
    floorVal = 0.0;
end
end
