function res = run_paper2_source_beam(cfg)
%RUN_PAPER2_SOURCE_BEAM Reproduce 1D source-beam transport benchmark (Section 6.1.2).

if nargin < 1 || isempty(cfg)
    cfg = mm_default_config();
end

if ~exist(cfg.paths.results, 'dir')
    mkdir(cfg.paths.results);
end

closeFigures = get_close_figures(cfg);
tCfg = get_timing_cfg(cfg);
parallel_cfg = get_parallel_cfg(cfg);
sb_cfg = get_source_beam_cfg(cfg);

tRef = tic;
ref_cfg = get_field_or(cfg, 'reference_source_beam', struct());
ref_cfg.parallel = parallel_cfg;
ref_cfg.domain = get_field_or(ref_cfg, 'domain', sb_cfg.domain);
ref_cfg.tf = get_field_or(ref_cfg, 'tf', sb_cfg.tf);
ref_cfg.psi_vac_density = get_field_or(ref_cfg, 'psi_vac_density', cfg.physics.psi_vac_density);
ref_cfg.beam_exponent = get_field_or(ref_cfg, 'beam_exponent', sb_cfg.beam_exponent);
ref_cfg.left_boundary_normalize = get_field_or(ref_cfg, 'left_boundary_normalize', sb_cfg.left_boundary_normalize);
ref_cfg.sigma_s = get_field_or(ref_cfg, 'sigma_s', get_field_or(sb_cfg, 'sigma_s', []));
ref_cfg.sigma_a = get_field_or(ref_cfg, 'sigma_a', get_field_or(sb_cfg, 'sigma_a', []));
ref_cfg.Q = get_field_or(ref_cfg, 'Q', get_field_or(sb_cfg, 'Q', []));
if isfield(sb_cfg, 'left_boundary_psi')
    ref_cfg.boundary_psi_left = sb_cfg.left_boundary_psi;
end
if isfield(sb_cfg, 'right_boundary_psi')
    ref_cfg.boundary_psi_right = sb_cfg.right_boundary_psi;
end
ref = mm_reference_sn_source_beam(ref_cfg);
timeRef = toc(tRef);

zL = sb_cfg.domain(1);
zR = sb_cfg.domain(2);
nCells = sb_cfg.n_cells;
edges = linspace(zL, zR, nCells + 1);
z = 0.5 * (edges(1:end-1) + edges(2:end));
dz = edges(2) - edges(1);
gridData = struct('z', z(:), 'dz', dz, 'nCells', nCells, 'ghost', 1);

rhoRef = interp1(ref.z, ref.rho, z, 'linear', 'extrap');

registry = mm_model_registry(cfg);
rows_cell = cell(numel(registry), 1);
run_cell = cell(numel(registry), 1);

ens_ctx = mm_parallel_context(parallel_cfg, numel(registry), 'ensemble');
use_ensemble = ens_ctx.use_parallel;

if use_ensemble
    parfor i = 1:numel(registry)
        rec = registry(i);
        out = run_single_source_beam_model(rec, cfg, sb_cfg, gridData, z, dz, nCells, rhoRef, timeRef, tCfg, parallel_cfg, false);
        rows_cell{i} = out.row;
        run_cell{i} = out;
    end
else
    for i = 1:numel(registry)
        rec = registry(i);
        out = run_single_source_beam_model(rec, cfg, sb_cfg, gridData, z, dz, nCells, rhoRef, timeRef, tCfg, parallel_cfg, true);
        rows_cell{i} = out.row;
        run_cell{i} = out;
        if tCfg.print_summary
            print_timing_summary(out.row);
        end
    end
end

if use_ensemble && tCfg.print_summary
    for i = 1:numel(run_cell)
        print_timing_summary(run_cell{i}.row);
    end
end

rows = vertcat(rows_cell{:});
T = struct2table(rows);
outCsv = fullfile(cfg.paths.results, 'paper2_source_beam_errors.csv');
writetable(T, outCsv);

plot_source_beam_errors(T, fullfile(cfg.paths.results, 'paper2_source_beam_errors'), closeFigures);

res = struct();
res.table = T;
res.csv = outCsv;
res.reference = ref;

end

function out = run_single_source_beam_model(rec, cfg, sb_cfg, gridData, z, dz, nCells, rhoRef, timeRef, tCfg, parallel_cfg, show_progress)
tModelBuild = tic;
model = mm_build_model(rec.name, rec.order, cfg.models);
timeModelBuild = toc(tModelBuild);

psi0 = cfg.physics.psi_vac_density * ones(1, nCells);
u = model.b_iso * psi0;

[sigma_s, sigma_a, Q] = source_beam_fields(z, sb_cfg);
beamExp = get_field_or(sb_cfg, 'beam_exponent', 1.0e5);
normalizeLeft = logical(get_field_or(sb_cfg, 'left_boundary_normalize', false));
left_profile = get_field_or(sb_cfg, 'left_boundary_psi', @(mu_vals, t) source_beam_left_profile(mu_vals, t, beamExp, normalizeLeft));
right_profile = get_field_or(sb_cfg, 'right_boundary_psi', @(mu_vals, t) cfg.physics.psi_vac_density * ones(size(mu_vals)));

phys = struct();
phys.sigma_s = sigma_s;
phys.sigma_a = sigma_a;
phys.Q = Q;
phys.psi_vac_density = cfg.physics.psi_vac_density;
phys.boundary = struct('left', cfg.physics.psi_vac_density, 'right', cfg.physics.psi_vac_density);
phys.boundary_psi = struct('left', left_profile, 'right', right_profile);

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
solver_cfg.flux.parallel = parallel_cfg;
solver_cfg.flux.quad_flux = mm_build_quadrature(model, cfg.quad, 'flux');
timeQuadBuild = toc(tQuadBuild);

cflSafety = sb_cfg.cfl_safety;
if model.is_partial && model.needs_entropy
    cflSafety = sb_cfg.cfl_safety_partial_entropy;
end
dtCFL = cflSafety * ((1 - cfg.optimizer.eps_gamma) / 2.0) * dz;
nSteps = ceil(sb_cfg.tf / dtCFL);
dt = sb_cfg.tf / nSteps;

state = struct();
tStart = tic;
sumStep = init_step_timing_acc();
for it = 1:nSteps
    [u, state] = mm_step_strang(u, model, phys, gridData, dt, solver_cfg, state);
    if show_progress
        fprintf('\r[%s-%d] step %d/%d (%.1f%%)', rec.name, rec.order, it, nSteps, 100 * it / nSteps);
        if it == nSteps
            fprintf('\n');
        end
    end
    if tCfg.enabled && tCfg.collect_step_breakdown
        sumStep = accumulate_step_timing(sumStep, state);
    end
end
runtime = toc(tStart);

rho = (model.alpha1.' * u).';

err = abs(rho - rhoRef(:));
L1 = trapz(z, err);
Linf = max(err);
minRho = min(rho);

row = struct();
row.model = rec.name;
row.order = rec.order;
row.nMom = rec.nMom;
row.L1 = L1;
row.Linf = Linf;
row.symmetry_L1_rel = NaN;
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

write_model_checkpoint(cfg.paths.results, rec, z, rho, rhoRef, row);
write_model_plot(cfg.paths.results, rec, z, rhoRef, rho);

out = struct();
out.rec = rec;
out.rho = rho;
out.row = row;
end

function write_model_checkpoint(results_dir, rec, z, rho, rhoRef, row)
outProf = fullfile(results_dir, sprintf('paper2_source_beam_%s_%d_profile.csv', rec.name, rec.order));
Tprof = table(z(:), rho(:), rhoRef(:), 'VariableNames', {'z', 'rho_model', 'rho_ref'});
writetable(Tprof, outProf);

outMetrics = fullfile(results_dir, sprintf('paper2_source_beam_%s_%d_metrics.csv', rec.name, rec.order));
writetable(struct2table(row), outMetrics);

outCheckpoint = fullfile(results_dir, sprintf('paper2_source_beam_%s_%d_checkpoint.mat', rec.name, rec.order));
save(outCheckpoint, 'row', 'rho');
end

function write_model_plot(results_dir, rec, z, rhoRef, rho)
fig = create_hidden_figure();
if isempty(fig)
    return;
end

plot(z, rhoRef, 'k-', 'LineWidth', 1.5); hold on;
plot(z, rho, 'r--', 'LineWidth', 1.2);
grid(gca, 'on');
xlabel('z'); ylabel('rho');
title(sprintf('Source beam: %s-%d', rec.name, rec.order));
legend('S_N reference', 'model', 'Location', 'best');
outBase = fullfile(results_dir, sprintf('paper2_source_beam_%s_%d', rec.name, rec.order));
write_figure_files(fig, outBase);
close(fig);
end

function print_timing_summary(row)
fprintf('[timing] %s-%d: total=%.3fs solver=%.3fs source=%.3fs flux=%.3fs jac=%.3fs rec=%.3fs lim+flux=%.3fs\n', ...
    row.model, row.order, row.runtime_s, row.time_step_sum_s, row.time_source_sum_s, ...
    row.time_flux_sum_s, row.time_flux_jacobian_sum_s, row.time_flux_reconstruct_sum_s, row.time_flux_limiter_flux_sum_s);
end

function plot_source_beam_errors(T, outBase, closeFigures)
modelOrder = {'HFMn', 'HFPn', 'PMMn', 'PMPn', 'MN', 'PN'};

presentMask = false(size(modelOrder));
for i = 1:numel(modelOrder)
    presentMask(i) = any(strcmp(string(T.model), string(modelOrder{i})));
end
models = modelOrder(presentMask);

extraModels = setdiff(cellstr(unique(string(T.model), 'stable')), models, 'stable');
models = [models, extraModels];

fig = create_hidden_figure();
if isempty(fig)
    return;
end

tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot_error_metric(ax1, T, models, 'L1', 'E_h^1', '(a) L^1 convergence');

ax2 = nexttile;
plot_error_metric(ax2, T, models, 'Linf', 'E_h^\infty', '(b) L^\infty convergence');

write_figure_files(fig, outBase);
if closeFigures
    close(fig);
else
    close(fig);
end
end

function plot_error_metric(ax, T, models, metricName, yLabelTxt, titleTxt)
hold(ax, 'on');
for i = 1:numel(models)
    m = models{i};
    tm = T(strcmp(string(T.model), string(m)), :);
    if isempty(tm)
        continue;
    end

    [x, idx] = sort(tm.nMom);
    y = tm.(metricName)(idx);
    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;
    x = x(valid);
    y = y(valid);
    if isempty(x)
        continue;
    end

    sty = model_curve_style(m);
    loglog(ax, x, y, 'Color', sty.color, 'LineStyle', sty.line, 'Marker', sty.marker, ...
        'LineWidth', 1.2, 'MarkerSize', 4.5, 'DisplayName', model_display_name(m));
end

grid(ax, 'on');
xlabel(ax, 'n');
ylabel(ax, yLabelTxt, 'Interpreter', 'tex');
title(ax, titleTxt, 'Interpreter', 'tex');
add_order_guides(ax);
legend(ax, 'Location', 'best');
end

function add_order_guides(ax)
xl = xlim(ax);
yl = ylim(ax);
if any(~isfinite([xl, yl])) || xl(1) <= 0 || yl(1) <= 0
    return;
end

xRef = logspace(log10(xl(1)), log10(xl(2)), 48);
x0 = sqrt(xl(1) * xl(2));
y0 = yl(1) * 3.0;
yOrd1 = y0 * (xRef / x0).^(-1);
yOrd2 = y0 * (xRef / x0).^(-2);

loglog(ax, xRef, yOrd1, '--', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 0.9, 'HandleVisibility', 'off');
loglog(ax, xRef, yOrd2, '-.', 'Color', [0.15, 0.15, 0.15], 'LineWidth', 0.9, 'HandleVisibility', 'off');

idx1 = max(2, round(0.72 * numel(xRef)));
idx2 = max(2, round(0.50 * numel(xRef)));
text(ax, xRef(idx1), yOrd1(idx1) * 1.10, '1st order', 'Color', [0.35, 0.35, 0.35], ...
    'FontSize', 9, 'Interpreter', 'tex');
text(ax, xRef(idx2), yOrd2(idx2) * 1.10, '2nd order', 'Color', [0.15, 0.15, 0.15], ...
    'FontSize', 9, 'Interpreter', 'tex');
end

function sty = model_curve_style(modelName)
name = char(modelName);
switch name
    case 'HFMn'
        sty = struct('color', [0.0, 0.4470, 0.7410], 'line', '-', 'marker', 's');
    case 'HFPn'
        sty = struct('color', [0.3010, 0.7450, 0.9330], 'line', '-', 'marker', '^');
    case 'PMMn'
        sty = struct('color', [0.9290, 0.6940, 0.1250], 'line', '-', 'marker', 'd');
    case 'PMPn'
        sty = struct('color', [0.45, 0.45, 0.45], 'line', '-', 'marker', 'x');
    case 'MN'
        sty = struct('color', [0.4660, 0.6740, 0.1880], 'line', '--', 'marker', 'o');
    case 'PN'
        sty = struct('color', [0.0, 0.5, 0.0], 'line', ':', 'marker', 'o');
    otherwise
        sty = struct('color', [0.1, 0.1, 0.1], 'line', '-', 'marker', 'o');
end
end

function name = model_display_name(modelName)
name = char(modelName);
switch name
    case 'MN'
        name = 'M_N';
    case 'PN'
        name = 'P_N';
    case 'HFMn'
        name = 'HFM_n';
    case 'HFPn'
        name = 'HFP_n';
    case 'PMMn'
        name = 'PMM_n';
    case 'PMPn'
        name = 'PMP_n';
end
end
function [sigma_s, sigma_a, Q] = source_beam_fields(z, sb_cfg)
sigma_a = eval_cell_field(get_field_or(sb_cfg, 'sigma_a', []), z, @default_sigma_a);
sigma_s = eval_cell_field(get_field_or(sb_cfg, 'sigma_s', []), z, @default_sigma_s);
Q = eval_cell_field(get_field_or(sb_cfg, 'Q', []), z, @default_Q);
end

function v = eval_cell_field(spec, z, default_fun)
if isempty(spec)
    v = default_fun(z);
elseif isscalar(spec)
    v = spec * ones(size(z));
elseif isa(spec, 'function_handle')
    v = spec(z);
else
    v = spec;
end
v = reshape(v, 1, []);
if numel(v) ~= numel(z)
    error('Source-beam field has invalid size.');
end
end

function sigma_a = default_sigma_a(z)
sigma_a = zeros(size(z));
sigma_a(z <= 2.0) = 1.0;
end

function sigma_s = default_sigma_s(z)
sigma_s = zeros(size(z));
sigma_s(z > 1.0 & z <= 2.0) = 2.0;
sigma_s(z > 2.0) = 10.0;
end

function q = default_Q(z)
q = zeros(size(z));
q(z >= 1.0 & z <= 1.5) = 0.5;
end

function psi = source_beam_left_profile(mu_vals, ~, beam_exp, normalize_flag)
psi = exp(-beam_exp * (mu_vals - 1.0).^2);
if normalize_flag
    muPos = mu_vals(mu_vals >= 0);
    psiPos = psi(mu_vals >= 0);
    if ~isempty(muPos)
        denom = trapz(muPos, psiPos);
        if denom > 0
            psi = psi / denom;
        end
    end
end
end

function sb_cfg = get_source_beam_cfg(cfg)
if isfield(cfg, 'source_beam')
    sb_cfg = cfg.source_beam;
else
    sb_cfg = struct();
end
sb_cfg.domain = get_field_or(sb_cfg, 'domain', [0.0, 3.0]);
sb_cfg.tf = get_field_or(sb_cfg, 'tf', 2.5);
sb_cfg.n_cells = get_field_or(sb_cfg, 'n_cells', 1000);
sb_cfg.beam_exponent = get_field_or(sb_cfg, 'beam_exponent', 1.0e5);
sb_cfg.left_boundary_normalize = logical(get_field_or(sb_cfg, 'left_boundary_normalize', false));
sb_cfg.cfl_safety = get_field_or(sb_cfg, 'cfl_safety', get_field_or(cfg.solver, 'cfl_safety', 0.9));
sb_cfg.cfl_safety_partial_entropy = get_field_or(sb_cfg, 'cfl_safety_partial_entropy', ...
    get_field_or(cfg.solver, 'cfl_safety_partial_entropy', sb_cfg.cfl_safety));
end

function fig = create_hidden_figure()
fig = [];
try
    fig = figure('Color', 'w', 'Visible', 'off');
catch
    try
        set(groot, 'defaultFigureVisible', 'off');
        fig = figure('Color', 'w');
        set(fig, 'Visible', 'off');
    catch
        fig = [];
    end
end
end

function write_figure_files(fig, outBase)
try
    saveas(fig, [outBase '.png']);
catch
    try
        print(fig, [outBase '.png'], '-dpng', '-r200');
    catch
    end
end
try
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
catch
    try
        print(fig, [outBase '.pdf'], '-dpdf', '-painters');
    catch
    end
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

function par_cfg = get_parallel_cfg(cfg)
if isfield(cfg, 'parallel')
    par_cfg = cfg.parallel;
else
    par_cfg = struct();
end
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
