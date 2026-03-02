function [state_next, out_state] = mm_step_holo_m1pn_strang(state_holo, models, phys, grid, dt, solver_cfg, state)
%MM_STEP_HOLO_M1PN_STRANG Strang step for coupled HOLO M1+PN system.
% Reference for the split ordering: Seminarquelle 2, Section 4.1 (Strang splitting).

if nargin < 7 || isempty(state)
    state = struct();
end
tStrang = tic;

opt_cfg = get_field_or(solver_cfg, 'optimizer', struct('eps_gamma', 1.0e-2));
eps_gamma = get_field_or(opt_cfg, 'eps_gamma', 1.0e-2);

dt_max = ((1 - eps_gamma) / 2.0) * grid.dz;
if dt > dt_max
    if get_field_or(solver_cfg, 'auto_clip_dt', true)
        dt = 0.99 * dt_max;
    else
        error('CFL-like condition violated: dt=%g > %g', dt, dt_max);
    end
end

model_pn = models.pn;
model_m1 = models.m1;

source_cfg = get_field_or(solver_cfg, 'source', struct());
flux_cfg = get_field_or(solver_cfg, 'flux', struct());

[u_pn_half, s1_pn] = mm_step_source_exact(state_holo.u_pn, model_pn, phys, 0.5 * dt, source_cfg);
[~, s1_m1] = mm_step_source_exact(state_holo.u_m1, model_m1, phys, 0.5 * dt, source_cfg);
u_m1_half = u_pn_half(1:2, :);

state_half = struct('u_pn', u_pn_half, 'u_m1', u_m1_half);
[state_flux, sf] = mm_step_holo_m1pn_flux_rk2(state_half, models, phys, grid, dt, flux_cfg, state);

[u_pn_next, s2_pn] = mm_step_source_exact(state_flux.u_pn, model_pn, phys, 0.5 * dt, source_cfg);
[~, s2_m1] = mm_step_source_exact(state_flux.u_m1, model_m1, phys, 0.5 * dt, source_cfg);
u_m1_next = u_pn_next(1:2, :);

state_next = struct();
state_next.u_pn = u_pn_next;
state_next.u_m1 = u_m1_next;
state_next.diag = get_field_or(state_flux, 'diag', struct());
state_next.cache_pn = get_field_or(state_flux, 'cache_pn', struct());
state_next.cache_m1 = get_field_or(state_flux, 'cache_m1', struct());

out_state = struct();
out_state.last_source_1 = struct('pn', s1_pn, 'm1', s1_m1);
out_state.last_flux = sf;
out_state.last_source_2 = struct('pn', s2_pn, 'm1', s2_m1);
out_state.last_dt = dt;
out_state.cache_pn = state_next.cache_pn;
out_state.cache_m1 = state_next.cache_m1;
out_state.timing = struct();
out_state.timing.total_s = toc(tStrang);
out_state.timing.source_1_pn_s = get_nested_timing(s1_pn);
out_state.timing.source_1_m1_s = get_nested_timing(s1_m1);
out_state.timing.source_2_pn_s = get_nested_timing(s2_pn);
out_state.timing.source_2_m1_s = get_nested_timing(s2_m1);
out_state.timing.flux_s = get_timing_field(sf, 'total_s');

end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function v = get_timing_field(s, name)
v = 0.0;
if isfield(s, 'timing') && isfield(s.timing, name)
    v = s.timing.(name);
end
end

function v = get_nested_timing(s)
v = get_timing_field(s, 'total_s');
end
