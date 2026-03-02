function res = run_paper2_plane_source(cfg)
%RUN_PAPER2_PLANE_SOURCE Reproduce 1D plane-source transport benchmark.
% Reference setup: Seminarquelle 2, Section 6.1.1 (plane-source test).

if nargin < 1 || isempty(cfg)
    cfg = mm_default_config();
end

if ~exist(cfg.paths.results, 'dir')
    mkdir(cfg.paths.results);
end

closeFigures = get_close_figures(cfg);
tCfg = get_timing_cfg(cfg);

tRef = tic;
ref = mm_reference_sn_plane_source(cfg.reference);
timeRef = toc(tRef);

zL = cfg.paper2.domain(1);
zR = cfg.paper2.domain(2);
nCells = cfg.paper2.n_cells;
edges = linspace(zL, zR, nCells + 1);
z = 0.5 * (edges(1:end-1) + edges(2:end));
dz = edges(2) - edges(1);
gridData = struct('z', z(:), 'dz', dz, 'nCells', nCells, 'ghost', 1);

rhoRef = interp1(ref.z, ref.rho, z, 'linear', 'extrap');

registry = mm_model_registry(cfg);
rows = [];

for i = 1:numel(registry)
    tModelBuild = tic;
    rec = registry(i);
    model = mm_build_model(rec.name, rec.order, cfg.models);
    timeModelBuild = toc(tModelBuild);

    % Plane-source initial data (Section 5.5 / 6.1.1):
    % isotropic vacuum plus split Dirac mass at z=0 across adjacent cells.
    psi0 = cfg.physics.psi_vac_density * ones(1, nCells);
    [~, iCenter] = min(abs(z));
    if z(iCenter) <= 0
        iL = iCenter;
        iR = min(nCells, iCenter + 1);
    else
        iR = iCenter;
        iL = max(1, iCenter - 1);
    end
    psi0(iL) = psi0(iL) + 1 / (2 * dz);
    psi0(iR) = psi0(iR) + 1 / (2 * dz);

    u = model.b_iso * psi0;

    phys = struct();
    phys.sigma_s = cfg.paper2.sigma_s;
    phys.sigma_a = cfg.paper2.sigma_a;
    phys.Q = cfg.paper2.Q;
    phys.psi_vac_density = cfg.physics.psi_vac_density;
    phys.boundary = struct('left', cfg.physics.psi_vac_density, 'right', cfg.physics.psi_vac_density);

    tQuadBuild = tic;
    solver_cfg = struct();
    solver_cfg.auto_clip_dt = true;
    solver_cfg.optimizer = cfg.optimizer;
    solver_cfg.source = struct();
    solver_cfg.flux = struct();
    solver_cfg.flux.optimizer = cfg.optimizer;
    solver_cfg.flux.reconstruction = cfg.reconstruction;
    solver_cfg.flux.limiter = cfg.limiter;
    solver_cfg.flux.quad_cfg = cfg.quad;
    solver_cfg.flux.quad_flux = mm_build_quadrature(model, cfg.quad, 'flux');
    timeQuadBuild = toc(tQuadBuild);

    % CFL-like condition used in the paper for realizability-preserving RK stages.
    dtCFL = cfg.solver.cfl_safety * ((1 - cfg.optimizer.eps_gamma) / 2.0) * dz;
    nSteps = ceil(cfg.paper2.tf / dtCFL);
    dt = cfg.paper2.tf / nSteps;

    state = struct();
    tStart = tic;
    sumStep = init_step_timing_acc();
    for it = 1:nSteps
        [u, state] = mm_step_strang(u, model, phys, gridData, dt, solver_cfg, state);
        if tCfg.enabled && tCfg.collect_step_breakdown
            sumStep = accumulate_step_timing(sumStep, state);
        end
    end
    runtime = toc(tStart);

    rho = (model.alpha1.' * u).';

    err = abs(rho - rhoRef(:));
    L1 = trapz(z, err);
    Linf = max(err);

    rhoFlip = flipud(rho);
    symErr = norm(rho - rhoFlip, 1) / max(norm(rho, 1), eps);

    minRho = min(rho);

    row = struct();
    row.model = rec.name;
    row.order = rec.order;
    row.nMom = rec.nMom;
    row.L1 = L1;
    row.Linf = Linf;
    row.symmetry_L1_rel = symErr;
    row.min_rho = minRho;
    row.runtime_s = runtime;
    row.n_steps = nSteps;
    row.time_reference_s = timeRef;
    row.time_model_build_s = timeModelBuild;
    row.time_quad_build_s = timeQuadBuild;
    row.time_solver_total_s = runtime;
    row.time_step_sum_s = sumStep.total_s;
    row.time_source_sum_s = sumStep.source_1_s + sumStep.source_2_s;
    row.time_flux_sum_s = sumStep.flux_s;
    row.time_flux_jacobian_sum_s = sumStep.flux_jacobian_s;
    row.time_flux_reconstruct_sum_s = sumStep.flux_reconstruct_s;
    row.time_flux_limiter_flux_sum_s = sumStep.flux_limiter_flux_s;
    rows = [rows; row]; %#ok<AGROW>
    if tCfg.print_summary
        fprintf('[timing] %s-%d: total=%.3fs solver=%.3fs source=%.3fs flux=%.3fs jac=%.3fs rec=%.3fs lim+flux=%.3fs\n', ...
            rec.name, rec.order, runtime, sumStep.total_s, sumStep.source_1_s + sumStep.source_2_s, ...
            sumStep.flux_s, sumStep.flux_jacobian_s, sumStep.flux_reconstruct_s, sumStep.flux_limiter_flux_s);
    end

    outProf = fullfile(cfg.paths.results, sprintf('paper2_plane_source_%s_%d_profile.csv', rec.name, rec.order));
    Tprof = table(z(:), rho(:), rhoRef(:), 'VariableNames', {'z', 'rho_model', 'rho_ref'});
    writetable(Tprof, outProf);

    fig = figure('Color', 'w');
    plot(z, rhoRef, 'k-', 'LineWidth', 1.5); hold on;
    plot(z, rho, 'r--', 'LineWidth', 1.2);
    grid(gca, 'on');
    xlabel('z'); ylabel('rho');
    title(sprintf('Plane source: %s-%d', rec.name, rec.order));
    legend('S_N reference', 'model', 'Location', 'best');
    outBase = fullfile(cfg.paths.results, sprintf('paper2_plane_source_%s_%d', rec.name, rec.order));
    saveas(fig, [outBase '.png']);
    try
        exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
    catch
    end
    if closeFigures
        close(fig);
    end
end

T = struct2table(rows);
outCsv = fullfile(cfg.paths.results, 'paper2_plane_source_errors.csv');
writetable(T, outCsv);

plot_plane_errors(T, fullfile(cfg.paths.results, 'paper2_plane_source_errors'), closeFigures);

res = struct();
res.table = T;
res.csv = outCsv;
res.reference = ref;

end

function plot_plane_errors(T, outBase, closeFigures)
models = unique(string(T.model), 'stable');

figure('Color', 'w');
tiledlayout(1, 2);

nexttile;
hold on;
for i = 1:numel(models)
    m = models(i);
    tm = T(strcmp(string(T.model), m), :);
    [x, idx] = sort(tm.nMom);
    y = tm.L1(idx);
    loglog(x, y, '-o', 'DisplayName', char(m), 'LineWidth', 1.3);
end
xlabel('nMom'); ylabel('L1 error'); title('Plane source L1 error'); legend('Location', 'best');
grid(gca, 'on');

nexttile;
hold on;
for i = 1:numel(models)
    m = models(i);
    tm = T(strcmp(string(T.model), m), :);
    [x, idx] = sort(tm.nMom);
    y = tm.Linf(idx);
    loglog(x, y, '-o', 'DisplayName', char(m), 'LineWidth', 1.3);
end
xlabel('nMom'); ylabel('Linf error'); title('Plane source Linf error'); legend('Location', 'best');
grid(gca, 'on');

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
    if isfield(cfg.timing, 'collect_step_breakdown')
        tCfg.collect_step_breakdown = logical(cfg.timing.collect_step_breakdown);
    end
    if isfield(cfg.timing, 'print_summary')
        tCfg.print_summary = logical(cfg.timing.print_summary);
    end
end
if ~tCfg.enabled
    tCfg.collect_step_breakdown = false;
    tCfg.print_summary = false;
end
end

function acc = init_step_timing_acc()
acc = struct();
acc.total_s = 0.0;
acc.source_1_s = 0.0;
acc.source_2_s = 0.0;
acc.flux_s = 0.0;
acc.flux_jacobian_s = 0.0;
acc.flux_reconstruct_s = 0.0;
acc.flux_limiter_flux_s = 0.0;
end

function acc = accumulate_step_timing(acc, state)
if ~isfield(state, 'timing')
    return;
end
acc.total_s = acc.total_s + get_field_or(state.timing, 'total_s', 0.0);
acc.source_1_s = acc.source_1_s + get_field_or(state.timing, 'source_1_s', 0.0);
acc.source_2_s = acc.source_2_s + get_field_or(state.timing, 'source_2_s', 0.0);
acc.flux_s = acc.flux_s + get_field_or(state.timing, 'flux_s', 0.0);
if isfield(state, 'last_flux') && isfield(state.last_flux, 'timing')
    ft = state.last_flux.timing;
    acc.flux_jacobian_s = acc.flux_jacobian_s + get_field_or(ft, 'jacobian_s', 0.0);
    acc.flux_reconstruct_s = acc.flux_reconstruct_s + get_field_or(ft, 'reconstruct_s', 0.0);
    acc.flux_limiter_flux_s = acc.flux_limiter_flux_s + get_field_or(ft, 'limiter_and_flux_s', 0.0);
end
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
