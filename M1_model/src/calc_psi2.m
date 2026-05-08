function psi2 = calc_psi2(u)
%CALC_PSI2 Second angular moment of the 1D M1 closure.
% The active M1 model stays closed through the Eddington factor
% chi = chi(|j| / rho). See docs/REFERENCES.md for the code-to-source
% mapping. The reduced flux is clamped to [0, 1) to keep the square root
% in the realizable M1 closure well-defined.

    % Avoid division by zero in near-vacuum cells.
    eps_rho = 1e-14;
    eps_f = 1e-12;

    % Reduced flux f = |j| / rho.
    f = abs(u(2,:))./(u(1,:) + eps_rho);
    f = min(max(f, 0), 1 - eps_f);

    % One-dimensional Eddington factor chi(f).
    eddington = @(fval) ((3 + 4 * fval.^2) ...
        ./ (5 + 2 * sqrt(max(4 - 3 * fval.^2, eps_f))));

    % Second moment m2 = chi(f) * rho.
    psi2 = eddington(f) .* u(1,:);
end
