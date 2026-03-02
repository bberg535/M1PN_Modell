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
% LP feasibility problem Bw=u, w>=0 (Eq. (5.26)).
B = quad.B.';
b = u;

if isfield(model, 'lp_options')
    options = model.lp_options;
else
    options = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
end

nQ = size(B, 2);
try
    [~, ~, exitflag] = linprog(zeros(nQ, 1), [], [], B, b, zeros(nQ, 1), [], options);
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
