function [uL_lim, uR_lim, lim_state] = mm_apply_realizability_limiter(uCell, uL, uR, model, lim_cfg, quad)
%MM_APPLY_REALIZABILITY_LIMITER Apply cell-wise limiter to reconstructed states.
% References (Seminarquelle 2):
% - Limited reconstruction u_i^theta: Eq. (4.17)
% - Minimal limiter search: Eq. (4.18)
% - LP-based component limiter in characteristic variables: Eq. (5.25)
% - LP feasibility check for realizability: Eq. (5.26)
% - PMM interval constraints with epsilon-shift: Eq. (5.28), based on Eq. (3.3)

[nMom, nCells] = size(uCell);

if nargin < 5 || isempty(lim_cfg)
    lim_cfg = struct();
end

epsR = get_field_or(lim_cfg, 'epsR', 0.0);
epsTheta = get_field_or(lim_cfg, 'eps_theta', 0.0);
maxBisect = get_field_or(lim_cfg, 'lp_bisect_iter', 32);
% NEW (optional): characteristic-component LP limiter for full-moment models.
useCharLP = logical(get_field_or(lim_cfg, 'paper_lp_characteristic', false));
par_cfg = get_field_or(lim_cfg, 'parallel', struct());
par_ctx = mm_parallel_context(par_cfg, nCells, 'cells');
use_par_cells = par_ctx.use_parallel;

uL_lim = uL;
uR_lim = uR;
thetaCell = zeros(1, nCells);
thetaComp = zeros(nMom, nCells);
usedCharLP = false(1, nCells);
fallbackCharLP = false(1, nCells);
charBasis = get_field_or(lim_cfg, 'characteristic_basis', struct());

if use_par_cells
    parfor i = 1:nCells
        [uL_i, uR_i, theta_i, thetaVec_i, usedChar_i, fallback_i] = ...
            limit_single_cell(i, uCell, uL, uR, model, quad, epsR, epsTheta, maxBisect, useCharLP, charBasis);
        uL_lim(:, i) = uL_i;
        uR_lim(:, i) = uR_i;
        thetaCell(i) = theta_i;
        thetaComp(:, i) = thetaVec_i;
        usedCharLP(i) = usedChar_i;
        fallbackCharLP(i) = fallback_i;
    end
else
    for i = 1:nCells
        [uL_i, uR_i, theta_i, thetaVec_i, usedChar_i, fallback_i] = ...
            limit_single_cell(i, uCell, uL, uR, model, quad, epsR, epsTheta, maxBisect, useCharLP, charBasis);
        uL_lim(:, i) = uL_i;
        uR_lim(:, i) = uR_i;
        thetaCell(i) = theta_i;
        thetaComp(:, i) = thetaVec_i;
        usedCharLP(i) = usedChar_i;
        fallbackCharLP(i) = fallback_i;
    end
end

lim_state = struct('theta', thetaCell);
lim_state.theta_components = thetaComp;
if useCharLP && strcmp(model.realizability, 'lp')
    lim_state.mode = 'paper_char_lp';
else
    lim_state.mode = 'paper_scalar';
end
lim_state.char_lp_used_cells = usedCharLP;
lim_state.char_lp_fallback_cells = fallbackCharLP;

end

function [uL_out, uR_out, theta, thetaVec, usedChar, fallbackChar] = limit_single_cell(i, uCell, uL, uR, model, quad, epsR, epsTheta, maxBisect, useCharLP, charBasis)
nMom = size(uCell, 1);
ubar = uCell(:, i);
uL_i = uL(:, i);
uR_i = uR(:, i);

usedChar = false;
fallbackChar = false;
theta = 1.0;
thetaVec = ones(nMom, 1);

[okBar, ~] = mm_is_realizable(ubar, model, quad, struct('epsR', epsR));
if okBar
    usedVectorLimiter = false;
    if useCharLP && strcmp(model.realizability, 'lp')
        [thetaVecTry, okVec] = theta_lp_characteristic_pair(ubar, uL_i, uR_i, model, quad, charBasis, i);
        if okVec
            thetaVecTry = min(1.0, max(0.0, thetaVecTry + epsTheta));
            [uL_try, uR_try] = apply_characteristic_theta(ubar, uL_i, uR_i, thetaVecTry, charBasis, i);

            [okL, ~] = mm_is_realizable(uL_try, model, quad, struct('epsR', epsR));
            [okR, ~] = mm_is_realizable(uR_try, model, quad, struct('epsR', epsR));
            if okL && okR
                uL_i = uL_try;
                uR_i = uR_try;
                thetaVec = thetaVecTry;
                theta = max(thetaVec);
                usedVectorLimiter = true;
                usedChar = true;
            else
                fallbackChar = true;
            end
        else
            fallbackChar = true;
        end
    end

    if ~usedVectorLimiter
        tL = theta_for_state(ubar, uL_i, model, epsR, quad, maxBisect);
        tR = theta_for_state(ubar, uR_i, model, epsR, quad, maxBisect);
        theta = max(tL, tR);
        thetaVec = theta * ones(nMom, 1);
    end
end

if ~usedChar
    theta = min(1.0, max(0.0, theta + epsTheta));
    uL_i = theta * ubar + (1 - theta) * uL_i;
    uR_i = theta * ubar + (1 - theta) * uR_i;
end

uL_out = uL_i;
uR_out = uR_i;
end

function theta = theta_for_state(ubar, urec, model, epsR, quad, maxBisect)
if strcmp(model.realizability, 'positive')
    theta = theta_positive(ubar, urec, epsR);
elseif strcmp(model.realizability, 'interval')
    theta = theta_interval(ubar, urec, model, epsR);
else
    theta = theta_lp(ubar, urec, model, quad, epsR, maxBisect);
end
end

function theta = theta_positive(ubar, urec, epsR)
theta = 0.0;
for j = 1:numel(ubar)
    if urec(j) < epsR
        den = ubar(j) - urec(j);
        if abs(den) > eps
            cand = (epsR - urec(j)) / den;
            if cand >= 0 && cand <= 1
                theta = max(theta, cand);
            else
                theta = 1.0;
            end
        else
            theta = 1.0;
        end
    end
end
theta = min(max(theta, 0.0), 1.0);
end

function theta = theta_interval(ubar, urec, model, epsR)
mu = model.mu_edges;
k = model.kIntervals;

if ~check_interval_constraints(ubar, mu, epsR)
    theta = 1.0;
    return;
end

theta = 0.0;
for j = 1:k
    a = mu(j);
    b = mu(j + 1);

    idx = 2*j - 1;
    u0b = ubar(idx);  u1b = ubar(idx + 1);
    u0r = urec(idx);  u1r = urec(idx + 1);

    sA = sqrt(a^2 + 1.0);
    sB = sqrt(b^2 + 1.0);

    theta = max(theta, theta_linear_bound(u0b, u0r, epsR));

    qb = u1b - a * u0b;
    qr = u1r - a * u0r;
    theta = max(theta, theta_linear_bound(qb, qr, epsR * sA));

    pb = b * u0b - u1b;
    pr = b * u0r - u1r;
    theta = max(theta, theta_linear_bound(pb, pr, epsR * sB));
end

theta = min(max(theta, 0.0), 1.0);
end

function ok = check_interval_constraints(u, mu, epsR)
k = (numel(u) / 2);
ok = true;
for j = 1:k
    a = mu(j);
    b = mu(j + 1);
    idx = 2*j - 1;
    u0 = u(idx);
    u1 = u(idx + 1);

    % Epsilon-shifted interval realizability constraints, Eq. (5.28).
    c1 = u0 - epsR;
    c2 = (u1 - a * u0) - epsR * sqrt(a^2 + 1.0);
    c3 = (b * u0 - u1) - epsR * sqrt(b^2 + 1.0);

    if min([c1, c2, c3]) <= 0
        ok = false;
        return;
    end
end
end

function theta = theta_linear_bound(vBar, vRec, bound)
if vRec >= bound
    theta = 0.0;
    return;
end

den = vBar - vRec;
if abs(den) <= eps
    theta = 1.0;
    return;
end
cand = (bound - vRec) / den;
if cand >= 0 && cand <= 1
    theta = cand;
else
    theta = 1.0;
end
end

function theta = theta_lp(ubar, urec, model, quad, epsR, maxBisect)
[okRec, ~] = mm_is_realizable(urec, model, quad, struct('epsR', epsR));
if okRec
    theta = 0.0;
    return;
end

lo = 0.0;
hi = 1.0;
for it = 1:maxBisect
    mid = 0.5 * (lo + hi);
    um = mid * ubar + (1 - mid) * urec;
    [okMid, ~] = mm_is_realizable(um, model, quad, struct('epsR', epsR));
    if okMid
        hi = mid;
    else
        lo = mid;
    end
end

theta = hi;
end

function [theta, ok] = theta_lp_characteristic_pair(ubar, uL, uR, model, quad, charBasis, iCell)
nMom = numel(ubar);

[V, Vinv] = get_characteristic_basis(charBasis, iCell, nMom);
if isempty(V)
    theta = ones(nMom, 1);
    ok = false;
    return;
end

uL = uL(:);
uR = uR(:);
ubar = ubar(:);

deltaL = Vinv * (ubar - uL);
deltaR = Vinv * (ubar - uR);

M1 = V * diag(deltaL);
M2 = V * diag(deltaR);
B = quad.B.';
nQ = size(B, 2);

% Eq. (5.25)-style LP for component-wise limiting in characteristic coordinates.
Aeq = [B, zeros(nMom, nQ), -M1; ...
       zeros(nMom, nQ), B, -M2];
beq = [uL; uR];

c = [zeros(2 * nQ, 1); ones(nMom, 1)];
lb = [zeros(2 * nQ, 1); zeros(nMom, 1)];
ub = [inf(2 * nQ, 1); ones(nMom, 1)];

if isfield(model, 'lp_options')
    options = model.lp_options;
else
    options = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
end

try
    [x, ~, exitflag] = linprog(c, [], [], Aeq, beq, lb, ub, options);
catch
    exitflag = -1;
    x = [];
end

ok = exitflag > 0 && ~isempty(x);
if ok
    theta = x((2 * nQ + 1):end);
else
    theta = ones(nMom, 1);
end
end

function [uL_lim, uR_lim] = apply_characteristic_theta(ubar, uL, uR, theta, charBasis, iCell)
nMom = numel(ubar);
[V, Vinv] = get_characteristic_basis(charBasis, iCell, nMom);
if isempty(V)
    uL_lim = uL;
    uR_lim = uR;
    return;
end

ubar_c = Vinv * ubar;
uL_c = Vinv * uL;
uR_c = Vinv * uR;

uL_lim = real(V * (theta .* ubar_c + (1 - theta) .* uL_c));
uR_lim = real(V * (theta .* ubar_c + (1 - theta) .* uR_c));
end

function [V, Vinv] = get_characteristic_basis(charBasis, iCell, nMom)
V = eye(nMom);
Vinv = eye(nMom);

if ~isstruct(charBasis)
    return;
end
if ~isfield(charBasis, 'V') || ~isfield(charBasis, 'Vinv')
    return;
end
if numel(charBasis.V) < iCell || numel(charBasis.Vinv) < iCell
    return;
end
if isempty(charBasis.V{iCell}) || isempty(charBasis.Vinv{iCell})
    return;
end

Vtry = charBasis.V{iCell};
VinvTry = charBasis.Vinv{iCell};

if any(~isfinite(Vtry(:))) || any(~isfinite(VinvTry(:))) || size(Vtry, 1) ~= nMom || size(Vtry, 2) ~= nMom
    return;
end

if ~isreal(Vtry) || ~isreal(VinvTry)
    if max(abs(imag(Vtry(:)))) > 1e-10 || max(abs(imag(VinvTry(:)))) > 1e-10
        return;
    end
    Vtry = real(Vtry);
    VinvTry = real(VinvTry);
end

V = Vtry;
Vinv = VinvTry;
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
