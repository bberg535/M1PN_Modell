function ref = sn_reference_periodic_source(cfg_ref)
%SN_REFERENCE_PERIODIC_SOURCE High-resolution S_N reference for the periodic source benchmark.
%   The setup matches the M1/M1PN scripts in this workspace:
%   - slab domain with periodic advection,
%   - isotropic scattering,
%   - isotropic cell source q/2,
%   - isotropic initial density rho/2.

    if nargin < 1
        cfg_ref = struct();
    end

    domain = get_field_or(cfg_ref, 'domain', [0.0, 2.0]);
    zL = domain(1);
    zR = domain(2);
    tf = get_field_or(cfg_ref, 'tf', 1.0);
    nCells = get_field_or(cfg_ref, 'n_cells', 1000);
    nMu = get_field_or(cfg_ref, 'n_mu', 128);
    cfl = get_field_or(cfg_ref, 'cfl', 0.45);

    edges = linspace(zL, zR, nCells + 1);
    z = 0.5 * (edges(1:end-1) + edges(2:end));
    dz = edges(2) - edges(1);

    sigma_a = eval_cell_field(get_field_or(cfg_ref, 'sigma_a', 0.0), z);
    sigma_s = eval_cell_field(get_field_or(cfg_ref, 'sigma_s', 0.0), z);
    q = eval_source_field(cfg_ref, z);
    rho0 = eval_initial_rho(cfg_ref, z, dz);

    [mu, w] = gausslegendre(nMu);
    max_speed = max(abs(mu));

    psi = 0.5 * (ones(nMu, 1) * rho0);

    t = 0.0;
    while t < tf - 1e-14
        dt = min(cfl * dz / max_speed, tf - t);

        k1 = rhs_sn(psi);
        psi1 = max(psi + dt * k1, 0.0);

        k2 = rhs_sn(psi1);
        psi = max(psi + 0.5 * dt * (k1 + k2), 0.0);

        t = t + dt;
    end

    rho = w.' * psi;
    j = (w(:) .* mu(:)).' * psi;
    m2 = (w(:) .* (mu(:) .^ 2)).' * psi;

    ref = struct();
    ref.z = z(:);
    ref.rho = rho(:);
    ref.j = j(:);
    ref.m2 = m2(:);
    ref.psi = psi;
    ref.q = q(:);
    ref.t = t;
    ref.mu = mu;
    ref.w = w;

    function rhs = rhs_sn(psi_state)
        rho_state = w.' * psi_state;
        src = bsxfun(@times, sigma_s, 0.5 * rho_state - psi_state) ...
            - bsxfun(@times, sigma_a, psi_state) ...
            + 0.5 * (ones(nMu, 1) * q);

        adv = zeros(size(psi_state));
        im1 = [nCells, 1:nCells-1];
        ip1 = [2:nCells, 1];
        for m = 1:nMu
            if mu(m) >= 0
                adv(m, :) = -mu(m) * (psi_state(m, :) - psi_state(m, im1)) / dz;
            else
                adv(m, :) = -mu(m) * (psi_state(m, ip1) - psi_state(m, :)) / dz;
            end
        end

        rhs = adv + src;
    end
end

function q = eval_source_field(cfg_ref, z)
    q_spec = get_field_or(cfg_ref, 'q', []);
    if isempty(q_spec)
        source_strength = get_field_or(cfg_ref, 'source_strength', 0.0);
        q_full = source_term(numel(z), source_strength);
        q = q_full(1, :);
    else
        q = eval_cell_field(q_spec, z);
    end
end

function rho0 = eval_initial_rho(cfg_ref, z, dz)
    rho_spec = get_field_or(cfg_ref, 'initial_rho', []);
    if isempty(rho_spec)
        rho0 = zeros(size(z));
        rho0(1) = 1.0 / dz;
    else
        rho0 = eval_cell_field(rho_spec, z);
    end
end

function v = eval_cell_field(spec, z)
    if isscalar(spec)
        v = spec * ones(size(z));
    elseif isa(spec, 'function_handle')
        v = spec(z);
    else
        v = spec;
    end

    v = reshape(v, 1, []);
    if numel(v) ~= numel(z)
        error('sn_reference_periodic_source:invalidFieldSize', ...
            'Expected field data of length %d.', numel(z));
    end
end

function v = get_field_or(s, name, default)
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end
