function res = run_holo_m1pn_plane_source(cfg)
%RUN_HOLO_M1PN_PLANE_SOURCE Plane-source benchmark for coupled HOLO M1+PN.

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

rows = [];
m1Family = get_field_or(cfg.holo, 'm1_family', 'MN');
m1Order = get_field_or(cfg.holo, 'm1_order', 1);

for i = 1:numel(cfg.holo.pn_orders)
    N = cfg.holo.pn_orders(i);
    tModelBuild = tic;
    model_pn = mm_build_model('PN', N, cfg.models);
    model_m1 = mm_build_model(m1Family, m1Order, cfg.models);
    timeModelBuild = toc(tModelBuild);

    if model_m1.nMom ~= 2
        error('Configured M1 model must have exactly 2 moments.');
    end

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

    u_pn0 = model_pn.b_iso * psi0;
    state_holo = struct('u_pn', u_pn0, 'u_m1', u_pn0(1:2, :), ...
        'cache_pn', struct(), 'cache_m1', struct(), 'diag', struct());

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
    solver_cfg.flux.quad_flux_pn = mm_build_quadrature(model_pn, cfg.quad, 'flux');
    solver_cfg.flux.quad_flux_m1 = mm_build_quadrature(model_m1, cfg.quad, 'flux');
    solver_cfg.flux.quad_lp_m1 = mm_build_quadrature(model_m1, cfg.quad, 'lp');
    solver_cfg.flux.holo = cfg.holo;
    solver_cfg.flux.sync = struct();
    solver_cfg.flux.sync.weight = get_field_or(cfg.holo, 'sync_weight', 'identity');
    solver_cfg.flux.sync.each_stage = logical(get_field_or(cfg.holo, 'sync_each_stage', true));
    timeQuadBuild = toc(tQuadBuild);

    dtCFL = cfg.solver.cfl_safety * ((1 - cfg.optimizer.eps_gamma) / 2.0) * dz;
    nSteps = ceil(cfg.paper2.tf / dtCFL);
    dt = cfg.paper2.tf / nSteps;

    iter_state = struct();
    thetaPnActive = 0;
    thetaM1Active = 0;
    syncStageSum = 0.0;
    syncFinalSum = 0.0;
    sumStep = init_holo_step_timing_acc();

    tStart = tic;
    for it = 1:nSteps
        [state_holo, iter_state] = mm_step_holo_m1pn_strang( ...
            state_holo, struct('pn', model_pn, 'm1', model_m1), phys, gridData, dt, solver_cfg, iter_state);

        if isfield(iter_state, 'last_flux') && isfield(iter_state.last_flux, 'diag')
            d = iter_state.last_flux.diag;
            thetaPnActive = thetaPnActive + get_field_or(d, 'theta_pn_active_stage1', 0) + get_field_or(d, 'theta_pn_active_stage2', 0);
            thetaM1Active = thetaM1Active + get_field_or(d, 'theta_m1_active_stage1', 0) + get_field_or(d, 'theta_m1_active_stage2', 0);
            syncStageSum = syncStageSum + get_field_or(d, 'sync_corr_stage1', 0.0);
            syncFinalSum = syncFinalSum + get_field_or(d, 'sync_corr_final', 0.0);
        end
        if tCfg.enabled && tCfg.collect_step_breakdown
            sumStep = accumulate_holo_step_timing(sumStep, iter_state);
        end
    end
    runtime = toc(tStart);

    rho = (model_pn.alpha1.' * state_holo.u_pn).';

    err = abs(rho - rhoRef(:));
    L1 = trapz(z, err);
    Linf = max(err);
    rhoFlip = flipud(rho);
    symErr = norm(rho - rhoFlip, 1) / max(norm(rho, 1), eps);
    minRho = min(rho);

    row = struct();
    row.model = 'HOLOM1PN';
    row.order = N;
    row.nMom = model_pn.nMom;
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
    row.time_flux_sum_s = sumStep.flux_s;
    row.time_source_pn_sum_s = sumStep.source_pn_s;
    row.time_source_m1_sum_s = sumStep.source_m1_s;
    row.time_flux_pn_limiter_flux_sum_s = sumStep.pn_limiter_flux_s;
    row.time_flux_m1_limiter_flux_sum_s = sumStep.m1_limiter_flux_s;
    row.time_flux_sync_sum_s = sumStep.sync_s;
    row.theta_pn_active_total = thetaPnActive;
    row.theta_m1_active_total = thetaM1Active;
    row.sync_corr_stage1_mean = syncStageSum / max(nSteps, 1);
    row.sync_corr_final_mean = syncFinalSum / max(nSteps, 1);
    rows = [rows; row]; %#ok<AGROW>
    if tCfg.print_summary
        fprintf('[timing] HOLOM1PN-%d: total=%.3fs solver=%.3fs source_pn=%.3fs source_m1=%.3fs flux=%.3fs pn_lim+flux=%.3fs m1_lim+flux=%.3fs sync=%.3fs\n', ...
            N, runtime, sumStep.total_s, sumStep.source_pn_s, sumStep.source_m1_s, sumStep.flux_s, ...
            sumStep.pn_limiter_flux_s, sumStep.m1_limiter_flux_s, sumStep.sync_s);
    end

    outProf = fullfile(cfg.paths.results, sprintf('paper2_plane_source_HOLOM1PN_%d_profile.csv', N));
    Tprof = table(z(:), rho(:), state_holo.u_m1(1, :).', rhoRef(:), ...
        'VariableNames', {'z', 'rho_model', 'rho_m1', 'rho_ref'});
    writetable(Tprof, outProf);

    fig = figure('Color', 'w');
    plot(z, rhoRef, 'k-', 'LineWidth', 1.5); hold on;
    plot(z, rho, 'r--', 'LineWidth', 1.2);
    grid(gca, 'on');
    xlabel('z'); ylabel('rho');
    title(sprintf('Plane source: HOLOM1PN-%d', N));
    legend('S_N reference', 'HOLOM1PN', 'Location', 'best');
    outBase = fullfile(cfg.paths.results, sprintf('paper2_plane_source_HOLOM1PN_%d', N));
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
outCsv = fullfile(cfg.paths.results, 'holo_m1pn_plane_source_errors.csv');
writetable(T, outCsv);

plot_holo_errors(T, fullfile(cfg.paths.results, 'holo_m1pn_plane_source_errors'), closeFigures);

res = struct();
res.table = T;
res.csv = outCsv;
res.reference = ref;

end

function plot_holo_errors(T, outBase, closeFigures)
figure('Color', 'w');
tiledlayout(1, 2);

nexttile;
[x1, idx1] = sort(T.nMom);
loglog(x1, T.L1(idx1), '-o', 'LineWidth', 1.3);
xlabel('nMom'); ylabel('L1 error'); title('HOLOM1PN plane source L1');
grid(gca, 'on');

nexttile;
[x2, idx2] = sort(T.nMom);
loglog(x2, T.Linf(idx2), '-o', 'LineWidth', 1.3);
xlabel('nMom'); ylabel('Linf error'); title('HOLOM1PN plane source Linf');
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

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
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

function acc = init_holo_step_timing_acc()
acc = struct();
acc.total_s = 0.0;
acc.flux_s = 0.0;
acc.source_pn_s = 0.0;
acc.source_m1_s = 0.0;
acc.pn_limiter_flux_s = 0.0;
acc.m1_limiter_flux_s = 0.0;
acc.sync_s = 0.0;
end

function acc = accumulate_holo_step_timing(acc, iterState)
if isfield(iterState, 'timing')
    t = iterState.timing;
    acc.total_s = acc.total_s + get_field_or(t, 'total_s', 0.0);
    acc.flux_s = acc.flux_s + get_field_or(t, 'flux_s', 0.0);
    acc.source_pn_s = acc.source_pn_s + get_field_or(t, 'source_1_pn_s', 0.0) + get_field_or(t, 'source_2_pn_s', 0.0);
    acc.source_m1_s = acc.source_m1_s + get_field_or(t, 'source_1_m1_s', 0.0) + get_field_or(t, 'source_2_m1_s', 0.0);
end
if isfield(iterState, 'last_flux') && isfield(iterState.last_flux, 'timing')
    tf = iterState.last_flux.timing;
    acc.pn_limiter_flux_s = acc.pn_limiter_flux_s + get_field_or(tf, 'pn_limiter_and_flux_s', 0.0);
    acc.m1_limiter_flux_s = acc.m1_limiter_flux_s + get_field_or(tf, 'm1_limiter_and_flux_s', 0.0);
    acc.sync_s = acc.sync_s + get_field_or(tf, 'sync_total_s', 0.0);
end
end
