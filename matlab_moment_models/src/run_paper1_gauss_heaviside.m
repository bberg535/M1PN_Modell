function res = run_paper1_gauss_heaviside(cfg)
%RUN_PAPER1_GAUSS_HEAVISIDE Reproduce slab ansatz tests from paper 1.

if nargin < 1 || isempty(cfg)
    cfg = mm_default_config();
end

if ~exist(cfg.paths.results, 'dir')
    mkdir(cfg.paths.results);
end

closeFigures = get_close_figures(cfg);
tCfg = get_timing_cfg(cfg);

muEval = cfg.paper1.eval_mu(:);
psiVac = cfg.physics.psi_vac_density;

psiGauss = @(mu) (1.0 / sqrt(2.0 * pi * cfg.paper1.gauss_sigma^2)) .* ...
    exp(-((mu - cfg.paper1.gauss_mu_bar).^2) ./ (2.0 * cfg.paper1.gauss_sigma^2));
psiHeav = @(mu) (mu < cfg.paper1.heaviside_jump) .* psiVac + (mu >= cfg.paper1.heaviside_jump) .* 1.0;

refGauss = psiGauss(muEval);
refHeav = psiHeav(muEval);

registry = mm_model_registry(cfg);
rows = [];

for i = 1:numel(registry)
    rec = registry(i);
    tModelBuild = tic;
    model = mm_build_model(rec.name, rec.order, cfg.models);
    timeModelBuild = toc(tModelBuild);
    tQuadBuild = tic;
    quad = mm_build_quadrature(model, cfg.quad, 'moment');
    timeQuadBuild = toc(tQuadBuild);

    [rowG, psiHatG] = evaluate_case('Gauss', psiGauss, refGauss, muEval, model, quad, cfg.optimizer, rec);
    [rowH, psiHatH] = evaluate_case('Heaviside', psiHeav, refHeav, muEval, model, quad, cfg.optimizer, rec);
    rowG.time_model_build_s = timeModelBuild;
    rowH.time_model_build_s = timeModelBuild;
    rowG.time_quad_build_s = timeQuadBuild;
    rowH.time_quad_build_s = timeQuadBuild;
    if tCfg.print_summary
        fprintf('[timing] %s-%d %s: total=%.3fs proj=%.3fs solve=%.3fs eval=%.3fs\n', rec.name, rec.order, rowG.case, rowG.runtime_s, rowG.time_project_s, rowG.time_solve_s, rowG.time_eval_s);
        fprintf('[timing] %s-%d %s: total=%.3fs proj=%.3fs solve=%.3fs eval=%.3fs\n', rec.name, rec.order, rowH.case, rowH.runtime_s, rowH.time_project_s, rowH.time_solve_s, rowH.time_eval_s);
    end

    rows = [rows; rowG; rowH]; %#ok<AGROW>

    outCase = fullfile(cfg.paths.results, sprintf('paper1_%s_%d_profiles.csv', rec.name, rec.order));
    profileTable = table(muEval, refGauss, psiHatG, refHeav, psiHatH, ...
        'VariableNames', {'mu', 'psi_ref_gauss', 'psi_hat_gauss', 'psi_ref_heaviside', 'psi_hat_heaviside'});
    writetable(profileTable, outCase);
end

T = struct2table(rows);
outCsv = fullfile(cfg.paths.results, 'paper1_gauss_heaviside_errors.csv');
writetable(T, outCsv);

plot_errors(T, 'Gauss', fullfile(cfg.paths.results, 'paper1_gauss_errors'), closeFigures);
plot_errors(T, 'Heaviside', fullfile(cfg.paths.results, 'paper1_heaviside_errors'), closeFigures);

res = struct();
res.table = T;
res.csv = outCsv;

end

function [row, psiHat] = evaluate_case(caseName, psiFun, ref, muEval, model, quad, opt_cfg, rec)
tStart = tic;

tProject = tic;
u = mm_project_density_to_moments(psiFun, model, quad);
tProjectTime = toc(tProject);
tEval = tic;
B_eval = model.basis_eval(muEval);
tEvalBasis = toc(tEval);

if model.needs_entropy
    tSolve = tic;
    [alpha, info] = mm_entropy_dual_solve(u, model, quad, opt_cfg, struct());
    tSolveTime = toc(tSolve);
    tEvalAnsatz = tic;
    psiHat = exp(min(B_eval * alpha, 700));
    tEvalAnsatzTime = toc(tEvalAnsatz);
    nIter = info.iterations;
    converged = info.converged;
else
    tSolve = tic;
    alpha = model.mass_matrix \ u;
    tSolveTime = toc(tSolve);
    tEvalAnsatz = tic;
    psiHat = B_eval * alpha;
    tEvalAnsatzTime = toc(tEvalAnsatz);
    nIter = 0;
    converged = true;
end

err = abs(psiHat - ref);
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
row.time_project_s = tProjectTime;
row.time_solve_s = tSolveTime;
row.time_eval_s = tEvalBasis + tEvalAnsatzTime;
row.iterations = nIter;
row.converged = converged;

end

function plot_errors(T, caseName, outBase, closeFigures)
mask = strcmp(T.case, caseName);
Tc = T(mask, :);
models = unique(string(Tc.model), 'stable');

figure('Color', 'w');
tiledlayout(1, 2);

nexttile;
hold on;
for i = 1:numel(models)
    m = models(i);
    tm = Tc(strcmp(string(Tc.model), m), :);
    [x, idx] = sort(tm.nMom);
    y = tm.L1(idx);
    loglog(x, y, '-o', 'DisplayName', char(m), 'LineWidth', 1.3);
end
xlabel('nMom'); ylabel('L1 error'); title(sprintf('%s: L1', caseName)); grid on; legend('Location', 'best');

nexttile;
hold on;
for i = 1:numel(models)
    m = models(i);
    tm = Tc(strcmp(string(Tc.model), m), :);
    [x, idx] = sort(tm.nMom);
    y = tm.Linf(idx);
    loglog(x, y, '-o', 'DisplayName', char(m), 'LineWidth', 1.3);
end
xlabel('nMom'); ylabel('Linf error'); title(sprintf('%s: Linf', caseName)); grid on; legend('Location', 'best');

saveas(gcf, [outBase '.png']);
try
    exportgraphics(gcf, [outBase '.pdf'], 'ContentType', 'vector');
catch
end
if closeFigures
    close(gcf);
end
end

function closeFigures = get_close_figures(cfg)
closeFigures = true;
if isfield(cfg, 'io') && isfield(cfg.io, 'close_figures')
    closeFigures = logical(cfg.io.close_figures);
end
end

function tCfg = get_timing_cfg(cfg)
tCfg = struct('enabled', true, 'collect_step_breakdown', true, 'print_summary', true);
if isfield(cfg, 'timing')
    if isfield(cfg.timing, 'enabled')
        tCfg.enabled = logical(cfg.timing.enabled);
    end
    if isfield(cfg.timing, 'print_summary')
        tCfg.print_summary = logical(cfg.timing.print_summary);
    end
end
if ~tCfg.enabled
    tCfg.print_summary = false;
end
end
