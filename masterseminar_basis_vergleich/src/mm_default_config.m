function cfg = mm_default_config()
%MM_DEFAULT_CONFIG Default configuration for slab-moment MATLAB project.

cfg = struct();

cfg.paths = struct();
cfg.paths.project_root = fileparts(fileparts(mfilename('fullpath')));
cfg.paths.results = fullfile(cfg.paths.project_root, 'results');

cfg.models = struct();
cfg.models.families = {'PN', 'MN', 'HFPn', 'HFMn', 'PMPn', 'PMMn'};
cfg.models.orders = struct();
cfg.models.orders.PN = [3, 7];
cfg.models.orders.MN = [3, 7];
cfg.models.orders.HFPn = [4, 8];
cfg.models.orders.HFMn = [4, 8];
cfg.models.orders.PMPn = [4, 8];
cfg.models.orders.PMMn = [4, 8];

cfg.quad = struct();
% Seminarquelle 2, Section 5.4:
% - first-order models: Gauss-Lobatto order 15
% - full moments in slab flux integrals: approx. order 2N+40 on half intervals
cfg.quad.lobatto_order = 15;
cfg.quad.base_legendre_order = 80;
cfg.quad.eval_points = 801;

cfg.optimizer = struct();
cfg.optimizer.tau = 1.0e-9;
cfg.optimizer.eps_gamma = 1.0e-2;
cfg.optimizer.max_iter = 200;
cfg.optimizer.armijo_c = 1.0e-3;
cfg.optimizer.armijo_backtrack = 0.5;
cfg.optimizer.step_min = 1.0e-12;
cfg.optimizer.grad_tol = 1.0e-11;
cfg.optimizer.regularization_r = [0, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 5.0e-2, 0.1, 0.5, 1.0];
cfg.optimizer.rho_vac = 1.0e-8;
cfg.optimizer.use_change_of_basis = true;
% PMMn tends to be more robust without Jacobian basis scaling in the dual Newton solve.
cfg.optimizer.use_change_of_basis_partial_entropy = false;
cfg.optimizer.change_of_basis_cond = 1.0e10;

cfg.reconstruction = struct();
cfg.reconstruction.use_characteristic = true;
% PMMn uses the Section-5.2 block generalized eigensolver for characteristic reconstruction.
cfg.reconstruction.use_characteristic_partial_entropy = true;

cfg.limiter = struct();
cfg.limiter.type = 'paper';
cfg.limiter.epsR = 1.0e-11;
cfg.limiter.eps_theta = 1.0e-11;
cfg.limiter.lp_bisect_iter = 32;
cfg.limiter.mcl_bisect_iter = 40;
cfg.limiter.mcl_lambda = 1.0;
% Optional paper-style LP limiter in characteristic components (Eq. 5.25-like).
cfg.limiter.paper_lp_characteristic = false;

cfg.solver = struct();
cfg.solver.cfl_safety = 0.9;
% Tighter CFL for PMMn-type entropy partial moments to reduce central overshoots.
cfg.solver.cfl_safety_partial_entropy = 0.9;

cfg.parallel = struct();
cfg.parallel.enabled = false;
cfg.parallel.mode = 'single_run';
cfg.parallel.num_workers = 0;
cfg.parallel.min_cells = 128;
cfg.parallel.min_interfaces = 128;
cfg.parallel.pool_type = 'local';
cfg.parallel.allow_fp_drift = true;
% PMMn entropy runs are sensitive; keep inner loops serial by default.
cfg.parallel.force_serial_partial_entropy = true;

cfg.physics = struct();
cfg.physics.psi_vac_density = 1.0e-8 / 2.0;

cfg.paper1 = struct();
cfg.paper1.gauss_sigma = 0.5;
cfg.paper1.gauss_mu_bar = 0.0;
cfg.paper1.heaviside_jump = 0.0;
cfg.paper1.crossing_a = 1.0e3; % Eq. 6.3
cfg.paper1.eval_mu = linspace(-1, 1, cfg.quad.eval_points);
cfg.paper1.figure3_use_cfg_registry = false;
cfg.paper1.figure3_orders = struct();
cfg.paper1.figure3_orders.PN = [1:2:49, 60];
cfg.paper1.figure3_orders.MN = [1:2:19, 29, 39, 49];
cfg.paper1.figure3_orders.HFPn = 2:2:60;
cfg.paper1.figure3_orders.HFMn = 2:2:60;
cfg.paper1.figure3_orders.PMPn = 2:2:60;
cfg.paper1.figure3_orders.PMMn = 2:2:60;
cfg.paper1.nonsmooth_reg_floor = 1.0e-2;
cfg.paper1.nonsmooth_reg_retry_floor = 5.0e-2;
cfg.paper1.nonsmooth_retry_spike_factor = 50.0;
cfg.paper1.plot_endpoint_trim = 1.0e-3;

cfg.paper2 = struct();
% Plane-source benchmark parameters (Seminarquelle 2, Section 6.1.1).
cfg.paper2.domain = [-1.2, 1.2];
cfg.paper2.tf = 1.0;
cfg.paper2.n_cells = 1000;
cfg.paper2.sigma_s = 1.0;
cfg.paper2.sigma_a = 0.0;
cfg.paper2.Q = 0.0;

cfg.source_beam = struct();
% Source-beam benchmark parameters (Seminarquelle 2, Section 6.1.2).
cfg.source_beam.domain = [0.0, 3.0];
cfg.source_beam.tf = 2.5;
cfg.source_beam.n_cells = 1000;
cfg.source_beam.beam_exponent = 1.0e5;
cfg.source_beam.left_boundary_normalize = false;
% Dedicated Source-Beam CFL (Section 6.1.2): choose ~1.73 to get about 971 steps for n_cells=1000.
cfg.source_beam.cfl_safety = 1.73;
cfg.source_beam.cfl_safety_partial_entropy = 1.73;

cfg.holo = struct();
cfg.holo.enabled = true;
cfg.holo.pn_orders = [3, 7];
cfg.holo.m1_family = 'MN';
cfg.holo.m1_order = 1;
cfg.holo.sync_each_stage = true;
cfg.holo.use_mcl_pn = true;
cfg.holo.use_mcl_m1 = true;
cfg.holo.sync_weight = 'identity';
cfg.holo.compare_include_baselines = true;

cfg.reference = struct();
cfg.reference.n_cells = 300;
cfg.reference.n_mu = 64;
cfg.reference.cfl = 0.45;
cfg.reference.sigma_s = cfg.paper2.sigma_s;
cfg.reference.sigma_a = cfg.paper2.sigma_a;
cfg.reference.Q = cfg.paper2.Q;
cfg.reference.domain = cfg.paper2.domain;
cfg.reference.tf = cfg.paper2.tf;
cfg.reference.psi_vac_density = cfg.physics.psi_vac_density;

cfg.reference_source_beam = struct();
cfg.reference_source_beam.n_cells = 1200;
cfg.reference_source_beam.n_mu = 128;
cfg.reference_source_beam.cfl = 0.45;
cfg.reference_source_beam.domain = cfg.source_beam.domain;
cfg.reference_source_beam.tf = cfg.source_beam.tf;
cfg.reference_source_beam.beam_exponent = cfg.source_beam.beam_exponent;
cfg.reference_source_beam.left_boundary_normalize = cfg.source_beam.left_boundary_normalize;
cfg.reference_source_beam.psi_vac_density = cfg.physics.psi_vac_density;

cfg.io = struct();
cfg.io.write_pdf = true;
cfg.io.close_figures = false;

cfg.timing = struct();
cfg.timing.enabled = true;
cfg.timing.collect_step_breakdown = true;
cfg.timing.print_summary = true;

end
