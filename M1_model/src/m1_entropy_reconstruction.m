function [psi, min_psi, beta] = m1_entropy_reconstruction(u, mu_eval)
%M1_ENTROPY_RECONSTRUCTION Reconstruct the 1D M1 minimum-entropy ansatz.
% Given the M1 moments (rho, j), reconstruct the continuous kinetic
% density psi(mu) = exp(alpha0 + alpha1 * mu) associated with the
% realizable M1 closure. The inversion uses the Langevin relation
% j/rho = coth(beta) - 1/beta with beta = alpha1.

    mu_eval = reshape(mu_eval, 1, []);
    rho = reshape(u(1, :), [], 1);
    j = reshape(u(2, :), [], 1);
    num_cells = numel(rho);

    psi = zeros(num_cells, numel(mu_eval));
    beta = zeros(num_cells, 1);
    eps_rho = 1.0e-14;
    eps_f = 1.0e-12;

    for idx = 1:num_cells
        rho_i = rho(idx);
        if rho_i <= eps_rho
            continue;
        end

        f_i = j(idx) / rho_i;
        f_i = min(max(f_i, -(1 - eps_f)), 1 - eps_f);
        beta_i = invert_langevin(f_i);
        beta(idx) = beta_i;
        psi(idx, :) = entropy_profile(rho_i, beta_i, mu_eval);
    end

    min_psi = min(psi, [], 2);
end

function beta = invert_langevin(f)
    if abs(f) < 1.0e-12
        beta = 0.0;
        return;
    end

    sign_f = sign(f);
    target = abs(f);
    left = 0.0;
    right = 1.0;
    while langevin(right) < target && right < 100.0
        right = 2.0 * right;
    end

    for iter = 1:80
        mid = 0.5 * (left + right);
        if langevin(mid) < target
            left = mid;
        else
            right = mid;
        end
    end

    beta = sign_f * 0.5 * (left + right);
end

function value = langevin(x)
    if abs(x) < 1.0e-10
        value = x / 3.0;
    else
        value = coth(x) - 1.0 / x;
    end
end

function psi_row = entropy_profile(rho, beta, mu_eval)
    if abs(beta) < 1.0e-10
        psi_row = 0.5 * rho * ones(size(mu_eval));
        return;
    end

    if beta > 0
        prefactor = rho * beta / (1.0 - exp(-2.0 * beta));
        psi_row = prefactor * exp(beta * (mu_eval - 1.0));
    else
        abs_beta = -beta;
        prefactor = rho * abs_beta / (1.0 - exp(-2.0 * abs_beta));
        psi_row = prefactor * exp(-abs_beta * (mu_eval + 1.0));
    end
end
