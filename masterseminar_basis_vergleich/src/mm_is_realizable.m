function [flag, margin] = mm_is_realizable(u, model, quad, lim_cfg)
%MM_IS_REALIZABLE Check realizability for a single moment vector.
% References:
% - HFM/PMM analytical conditions: Seminarquelle 2 Eq. (3.3), Eq. (5.28)
% - Full-moment numerical realizability via LP: Seminarquelle 2 Eq. (5.26)

if nargin < 4 || isempty(lim_cfg)
    lim_cfg = struct();
end
epsR = get_field_or(lim_cfg, 'epsR', 0.0);

u = u(:);

switch model.realizability
    case 'positive'
        minv = min(u);
        flag = minv > epsR;
        margin = minv - epsR;

    case 'interval'
        [flag, margin] = check_partial_interval(u, model, epsR);

    case 'lp'
        [flag, margin] = check_lp_realizable(u, model, quad);

    otherwise
        error('Unknown realizability type: %s', model.realizability);
end

end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function [flag, margin] = check_partial_interval(u, model, epsR)
mu = model.mu_edges;
k = model.kIntervals;
mins = inf(3*k, 1);
idx = 1;
for j = 1:k
    a = mu(j);
    b = mu(j + 1);
    u0 = u(2*j - 1);
    u1 = u(2*j);

    % Eq. (5.28): shifted constraints derived from Eq. (3.3).
    c1 = u0 - epsR;
    c2 = (u1 - a * u0) - epsR * sqrt(a^2 + 1.0);
    c3 = (b * u0 - u1) - epsR * sqrt(b^2 + 1.0);

    mins(idx) = c1; idx = idx + 1;
    mins(idx) = c2; idx = idx + 1;
    mins(idx) = c3; idx = idx + 1;
end
margin = min(mins);
flag = margin > 0;
end

function [flag, margin] = check_lp_realizable(u, model, quad)
if strcmp(get_field_or(model, 'family', ''), 'full') && isfield(quad, 'lp_hull')
    [flag, margin] = check_full_polytope(u, quad.lp_hull);
    return;
end

% LP feasibility problem Bw=u, w>=0 (Eq. (5.26)).
A = get_lp_matrix(quad);
b = u(:);

if isfield(model, 'lp_options')
    options = model.lp_options;
else
    options = linprog_options();
end

nQ = size(A, 2);
try
    [~, ~, exitflag] = linprog(zeros(nQ, 1), [], [], A, b, zeros(nQ, 1), [], options);
catch
    % If optimization toolbox is unavailable in runtime, mark non-realizable.
    exitflag = -1;
end

flag = exitflag > 0;
if flag
    margin = 1.0;
else
    margin = -1.0;
end
end

function [flag, margin] = check_full_polytope(u, hull)
u = u(:);
u0 = u(1);
tol = 1.0e-12 * max(1.0, norm(u, Inf));
if u0 <= tol
    flag = false;
    margin = u0 - tol;
    return;
end

if hull.dimension == 1
    x = u(2) / u0;
    margin = min([x - hull.bounds(1), hull.bounds(2) - x]);
    flag = margin >= -tol;
    return;
end

x = u(2:end) / u0;
viol = hull.normals * x - hull.offsets;
margin = -max(viol);
flag = margin >= -tol;
end

function A = get_lp_matrix(quad)
if isfield(quad, 'Bt_full') && ~isempty(quad.Bt_full)
    A = quad.Bt_full;
else
    A = full(quad.B.');
end
end

function options = linprog_options()
if exist('optimoptions', 'file') == 2
    options = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
else
    options = optimset('Display', 'off');
end
end
