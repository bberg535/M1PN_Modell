function cfg = m1pn_solver_config(opts)
%M1PN_SOLVER_CONFIG Normalize solver options for simulate_m1pn.

    if nargin < 1 || ~isstruct(opts)
        opts = struct();
    end

    cfg = struct();
    cfg.boundary = get_field_or(opts, 'boundary', 'periodic');
    cfg.psi_boundary = get_field_or(opts, 'psi_boundary', 0.0);
    cfg.pn_ho_method = get_field_or(opts, 'pn_ho_method', 1);
    cfg.coupling_mode = normalize_mode(get_field_or(opts, ...
        'coupling_mode', 'explicit_synced'));

    default_sync = strcmp(cfg.coupling_mode, 'explicit_synced') || ...
        strcmp(cfg.coupling_mode, 'implicit_holo');
    cfg.sync_auxiliary_pn = logical(get_field_or(opts, ...
        'sync_auxiliary_pn', default_sync));

    cfg.holo_outer_max_iters = get_positive_integer(opts, ...
        'holo_outer_max_iters', 12);
    cfg.holo_inner_max_iters = get_positive_integer(opts, ...
        'holo_inner_max_iters', 12);
    cfg.holo_tol = get_positive_scalar(opts, 'holo_tol', 1.0e-4);
    cfg.holo_relaxation = get_positive_scalar(opts, 'holo_relaxation', 1.0);
    if cfg.holo_relaxation > 1.0
        error('m1pn_solver_config:invalidRelaxation', ...
            'Expected holo_relaxation to be in (0, 1].');
    end

    cfg.holo_initial_guess = normalize_initial_guess(get_field_or(opts, ...
        'holo_initial_guess', 'explicit_predictor'));
end

function mode = normalize_mode(raw_mode)
    if isstring(raw_mode)
        raw_mode = char(raw_mode);
    end
    mode = lower(strtrim(raw_mode));

    switch mode
        case {'explicit_failsafe', 'current', 'default'}
            mode = 'explicit_failsafe';
        case {'explicit_synced', 'explicit_sync', 'synced'}
            mode = 'explicit_synced';
        case {'implicit_holo', 'holo_implicit', 'implicit'}
            mode = 'implicit_holo';
        otherwise
            error('m1pn_solver_config:unsupportedMode', ...
                'Unsupported coupling mode "%s".', raw_mode);
    end
end

function guess_mode = normalize_initial_guess(raw_guess)
    if isstring(raw_guess)
        raw_guess = char(raw_guess);
    end
    guess_mode = lower(strtrim(raw_guess));

    switch guess_mode
        case {'explicit_predictor', 'predictor'}
            guess_mode = 'explicit_predictor';
        case {'previous_state', 'previous'}
            guess_mode = 'previous_state';
        otherwise
            error('m1pn_solver_config:unsupportedInitialGuess', ...
                'Unsupported holo_initial_guess "%s".', raw_guess);
    end
end

function value = get_positive_integer(opts, field_name, default_value)
    value = get_field_or(opts, field_name, default_value);
    if ~isscalar(value) || value < 1 || abs(value - round(value)) > 0
        error('m1pn_solver_config:invalidInteger', ...
            'Expected %s to be a positive integer.', field_name);
    end
end

function value = get_positive_scalar(opts, field_name, default_value)
    value = get_field_or(opts, field_name, default_value);
    if ~isscalar(value) || ~isfinite(value) || value <= 0
        error('m1pn_solver_config:invalidScalar', ...
            'Expected %s to be a positive finite scalar.', field_name);
    end
end

function value = get_field_or(opts, field_name, default_value)
    if isfield(opts, field_name)
        value = opts.(field_name);
    else
        value = default_value;
    end
end
