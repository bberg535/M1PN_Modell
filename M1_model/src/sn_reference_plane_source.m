function ref = sn_reference_plane_source(cfg_ref)
%SN_REFERENCE_PLANE_SOURCE High-resolution S_N reference for the plane-source benchmark.
%   Benchmark parameters follow Section 6.1.1 of "Seminarquelle 2":
%   - domain [-1.2, 1.2],
%   - final time 1,
%   - isotropic split Dirac initial condition,
%   - vacuum boundary values,
%   - sigma_s = 1, sigma_a = 0, Q = 0.
%   The initialization convention can be switched via cfg_ref.init_mode:
%   - 'physical'     : rho has total unit mass (psi increment 1/(4*dz)),
%   - 'master_basis' : matches masterseminar_basis_vergleich exactly
%                      (psi increment 1/(2*dz)).

    if nargin < 1
        cfg_ref = struct();
    end

    domain = get_field_or(cfg_ref, 'domain', [-1.2, 1.2]);
    zL = domain(1);
    zR = domain(2);
    tf = get_field_or(cfg_ref, 'tf', 1.0);
    nCells = get_field_or(cfg_ref, 'n_cells', 2400);
    nMu = get_field_or(cfg_ref, 'n_mu', 256);
    cfl = get_field_or(cfg_ref, 'cfl', 0.45);
    psi_vac = get_field_or(cfg_ref, 'psi_vac_density', 0.5e-8);
    init_mode = lower(get_field_or(cfg_ref, 'init_mode', 'physical'));

    if mod(nCells, 2) ~= 0
        error('sn_reference_plane_source:oddCellCount', ...
            'Plane-source benchmark requires an even number of cells.');
    end

    edges = linspace(zL, zR, nCells + 1);
    z = 0.5 * (edges(1:end-1) + edges(2:end));
    dz = edges(2) - edges(1);

    sigma_a = eval_cell_field(get_field_or(cfg_ref, 'sigma_a', 0.0), z);
    sigma_s = eval_cell_field(get_field_or(cfg_ref, 'sigma_s', 1.0), z);
    q = eval_cell_field(get_field_or(cfg_ref, 'Q', 0.0), z);

    [mu, w] = gausslegendre(nMu);
    max_speed = max(abs(mu));

    psi = psi_vac * ones(nMu, nCells);
    [~, iCenter] = min(abs(z));
    if z(iCenter) <= 0
        iL = iCenter;
        iR = min(nCells, iCenter + 1);
    else
        iR = iCenter;
        iL = max(1, iCenter - 1);
    end
    psi_increment = split_dirac_increment(dz, init_mode);
    psi(:, iL) = psi(:, iL) + psi_increment;
    psi(:, iR) = psi(:, iR) + psi_increment;

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
    ref.t = t;
    ref.mu = mu;
    ref.w = w;
    ref.init_mode = init_mode;

    function rhs = rhs_sn(psi_state)
        rho_state = w.' * psi_state;
        src = bsxfun(@times, sigma_s, 0.5 * rho_state - psi_state) ...
            - bsxfun(@times, sigma_a, psi_state) ...
            + repmat(q, nMu, 1);

        adv = zeros(size(psi_state));
        pos = mu >= 0;
        neg = ~pos;

        if any(pos)
            left_state = [psi_vac * ones(sum(pos), 1), psi_state(pos, 1:end-1)];
            adv(pos, :) = -bsxfun(@times, mu(pos), (psi_state(pos, :) - left_state) / dz);
        end

        if any(neg)
            right_state = [psi_state(neg, 2:end), psi_vac * ones(sum(neg), 1)];
            adv(neg, :) = -bsxfun(@times, mu(neg), (right_state - psi_state(neg, :)) / dz);
        end

        rhs = adv + src;
    end
end

function psi_increment = split_dirac_increment(dz, init_mode)
    switch init_mode
        case 'physical'
            % rho is the angular integral over [-1, 1], so an isotropic
            % density increment rho_cell is represented by psi = rho_cell / 2.
            psi_increment = 1.0 / (4.0 * dz);

        case 'master_basis'
            % Matches masterseminar_basis_vergleich exactly.
            psi_increment = 1.0 / (2.0 * dz);

        otherwise
            error('sn_reference_plane_source:unknownInitMode', ...
                'Unknown plane-source init_mode "%s".', init_mode);
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
        error('sn_reference_plane_source:invalidFieldSize', ...
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
