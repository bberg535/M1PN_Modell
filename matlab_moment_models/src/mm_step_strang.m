function [u_next, state] = mm_step_strang(u, model, phys, grid, dt, solver_cfg, state)
%MM_STEP_STRANG One Strang-splitting step: source(dt/2)-flux(dt)-source(dt/2).
% Reference: Seminarquelle 2, Section 4.1 (Strang split of source and flux subsystems).

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

source_cfg = get_field_or(solver_cfg, 'source', struct());
flux_cfg = get_field_or(solver_cfg, 'flux', struct());

[u_half, s1] = mm_step_source_exact(u, model, phys, 0.5 * dt, source_cfg);
[u_flux, sf] = mm_step_flux_rk2(u_half, model, phys, grid, dt, flux_cfg, state);
[u_next, s2] = mm_step_source_exact(u_flux, model, phys, 0.5 * dt, source_cfg);

state.last_source_1 = s1;
state.last_flux = sf;
state.last_source_2 = s2;
state.last_dt = dt;
state.timing = struct();
state.timing.total_s = toc(tStrang);
state.timing.source_1_s = get_timing_field(s1, 'total_s');
state.timing.flux_s = get_timing_field(sf, 'total_s');
state.timing.source_2_s = get_timing_field(s2, 'total_s');

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
