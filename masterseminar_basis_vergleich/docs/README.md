# MATLAB Moment Models (Slab v1)

This subproject implements a slab-geometry reimplementation of the two Schneider/Leibner papers:

- Paper 1: ansatz reconstruction + realizability (Gauss/Heaviside)
- Paper 2: second-order realizability-preserving splitting idea (Strang + RK2 + limiter)

## Layout

- `src/`: core APIs and helpers
- `scripts/run_all.m`: run both benchmark scripts
- `tests/`: MATLAB unit tests (`runtests_mm`)
- `results/`: CSV and plot outputs

## Main API

- `cfg = mm_default_config()`
- `model = mm_build_model(model_name, order, cfg_model)`
- `quad = mm_build_quadrature(model, cfg_quad, purpose)`
- `u = mm_project_density_to_moments(psi_handle, model, quad)`
- `[alpha, info] = mm_entropy_dual_solve(u, model, quad, opt_cfg, cache_state)`
- `[V, Vinv, state] = mm_characteristic_basis(u, model, quad_flux, opt_cfg, cache_state)`
- `[uL, uR, rec_state] = mm_reconstruct_characteristic(u_cell, jac_state, grid, rec_cfg)`
- `[uL_lim, uR_lim, lim_state] = mm_apply_realizability_limiter(uCell, uL, uR, model, lim_cfg, quad)`
- `[flux, flux_state] = mm_kinetic_flux(uL, uR, model, quad, flux_cfg)`
- `[u_next, step_state] = mm_step_flux_rk2(u, model, phys, grid, dt, step_cfg, cache_state)`
- `[u_next, source_state] = mm_step_source_exact(u, model, phys, dt, source_cfg)`
- `[u_next, state] = mm_step_strang(u, model, phys, grid, dt, solver_cfg, state)`
- `[u_pn_sync, sync_state] = mm_sync_pn_to_m1(u_pn_ho, u_m1_target, sync_cfg)`
- `[state_next, step_state] = mm_step_holo_m1pn_flux_rk2(state, models, phys, grid, dt, step_cfg, cache_state)`
- `[state_next, out_state] = mm_step_holo_m1pn_strang(state, models, phys, grid, dt, solver_cfg, state)`
- `res = run_paper1_gauss_heaviside(cfg)`
- `res = run_paper1_figures_3_6(cfg)`
- `res = run_paper2_plane_source(cfg)`
- `res = run_holo_m1pn_plane_source(cfg)`
- `cmp = compare_holo_m1pn_vs_all(cfg)`
- `ref = mm_reference_sn_plane_source(cfg_ref)`

## Run

From MATLAB:

```matlab
addpath('matlab_moment_models/src');
addpath('matlab_moment_models/scripts');
run_all;
```

Limiter selection:

```matlab
cfg = mm_default_config();
cfg.limiter.type = 'paper'; % or 'mcl'
cfg.reconstruction.use_characteristic_partial_entropy = true; % PMMn via Section-5.2 block eigensolver
cfg.limiter.paper_lp_characteristic = false; % legacy scalar paper limiter (default)
% cfg.limiter.paper_lp_characteristic = true;  % NEW: characteristic-component LP limiter for PN/MN
cfg.io.close_figures = true; % false => figures stay open
res = run_paper2_plane_source(cfg);
```

Paper 1 (Figure 3-6 style slab tests):

```matlab
cfg = mm_default_config();
cfg.io.close_figures = false;  % keep figures open
cfg.paper1.figure3_use_cfg_registry = false; % use paper-like order sweep
res = run_paper1_figures_3_6(cfg);
```

Or via script:

```matlab
addpath('matlab_moment_models/scripts');
reproduce_paper1_figures_3_6
```

What is new (paper limiter path):

- Characteristic reconstruction now computes Section-5.2 eigenvectors from `J z = lambda H z` instead of finite-difference Jacobians; PMMn uses 2x2 block generalized eigenproblems.
- `paper_lp_characteristic = false`: old scalar limiter path (unchanged behavior)
- `paper_lp_characteristic = true`: NEW LP limiter in characteristic components for `PN/MN` (Eq.-5.25 style), with automatic fallback to the old scalar limiter if the LP is infeasible

Limiter comparison:

```matlab
addpath('matlab_moment_models/scripts');
compare_limiters('quick');      % fast comparison
compare_limiters('extended');   % broader comparison
compare_limiters('quick', true); % keep figures open

compare_holo('quick');          % HOLOM1PN vs baselines (quick)
compare_holo('extended');       % HOLOM1PN vs baselines (extended)
compare_holo('quick', true);    % keep figures open
```

Or tests:

```matlab
addpath('matlab_moment_models/tests');
runtests_mm;
```
