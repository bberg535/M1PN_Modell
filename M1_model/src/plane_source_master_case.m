function cfg = plane_source_master_case(varargin)
%PLANE_SOURCE_MASTER_CASE Plane-source setup matching masterseminar_basis_vergleich.
% The benchmark convention here intentionally follows the old
% masterseminar_basis_vergleich project:
% - domain [-1.2, 1.2],
% - n_cells = 1000,
% - psi_vac = 1e-8 / 2,
% - split Dirac initialization in the angular density psi,
% - reference grid 300 x 64 with interpolation back to the model grid.

    cfg = struct();
    cfg.domain = [-1.2, 1.2];
    cfg.T = 1.0;
    cfg.num_cells = 1000;
    cfg.psi_vac = 0.5e-8;
    cfg.sigma_a_value = 0.0;
    cfg.sigma_s_value = 1.0;
    cfg.source_strength = 0.0;

    % Matches mm_default_config: solver.cfl_safety and optimizer.eps_gamma.
    cfg.cfl_safety = 0.9;
    cfg.eps_gamma = 1.0e-2;

    % Matches the legacy S_N benchmark reference resolution.
    cfg.reference_n_cells = 300;
    cfg.reference_n_mu = 64;
    cfg.reference_cfl = 0.45;
    cfg.reference_init_mode = 'master_basis';

    cfg = apply_overrides(cfg, varargin{:});

    domain_length = cfg.domain(2) - cfg.domain(1);
    cfg.dz = domain_length / cfg.num_cells;
    cfg.z = linspace(cfg.domain(1) + 0.5 * cfg.dz, ...
        cfg.domain(2) - 0.5 * cfg.dz, cfg.num_cells).';

    dt_cfl = cfg.cfl_safety * ((1 - cfg.eps_gamma) / 2.0) * cfg.dz;
    cfg.num_steps = ceil(cfg.T / dt_cfl);
    cfg.dt = cfg.T / cfg.num_steps;

    cfg.rho0 = build_master_density(cfg.z, cfg.dz, cfg.psi_vac);
    cfg.sigma_a = cfg.sigma_a_value * ones(1, cfg.num_cells);
    cfg.sigma_s = cfg.sigma_s_value * ones(1, cfg.num_cells);
    cfg.vacuum_opts = struct('boundary', 'vacuum', 'psi_boundary', cfg.psi_vac);

    cfg.reference = struct();
    cfg.reference.domain = cfg.domain;
    cfg.reference.tf = cfg.T;
    cfg.reference.n_cells = cfg.reference_n_cells;
    cfg.reference.n_mu = cfg.reference_n_mu;
    cfg.reference.cfl = cfg.reference_cfl;
    cfg.reference.sigma_a = cfg.sigma_a_value;
    cfg.reference.sigma_s = cfg.sigma_s_value;
    cfg.reference.Q = cfg.source_strength;
    cfg.reference.psi_vac_density = cfg.psi_vac;
    cfg.reference.init_mode = cfg.reference_init_mode;
end

function rho0 = build_master_density(z, dz, psi_vac)
    num_cells = numel(z);
    rho0 = (2 * psi_vac) * ones(num_cells, 1);

    [~, i_center] = min(abs(z));
    if z(i_center) <= 0
        i_left = i_center;
        i_right = min(num_cells, i_center + 1);
    else
        i_right = i_center;
        i_left = max(1, i_center - 1);
    end

    % masterseminar_basis_vergleich initializes psi with 1/(2*dz) per cell.
    % Since rho = <psi> = 2 * psi for isotropic data, the corresponding rho
    % increment is 1/dz in each of the two center cells.
    rho0(i_left) = rho0(i_left) + 1.0 / dz;
    rho0(i_right) = rho0(i_right) + 1.0 / dz;
end

function cfg = apply_overrides(cfg, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('plane_source_master_case:invalidOverrides', ...
            'Overrides must be provided as name/value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k + 1};
        if ~isfield(cfg, name)
            error('plane_source_master_case:unknownField', ...
                'Unknown override field "%s".', name);
        end
        cfg.(name) = value;
    end
end
