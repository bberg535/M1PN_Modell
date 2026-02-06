function psi = scatter_halfstep(psi, w, tau, sigma_s, sigma_a, Q0)
% Exact update on each cell:
% psi' = -sigma_t psi + sigma_s*(phi/2) + Q0/2, with phi = sum w psi

sigma_t = sigma_s + sigma_a;

phi = (w.' * psi);                    % (1 x Nx)
psi_iso = (phi/2);                    % isotropic value per angle
psi_iso_mat = ones(size(psi,1),1) * psi_iso;

ea = exp(-sigma_a*tau);
es = exp(-sigma_s*tau);

psi = ea*( es*psi + (1-es)*psi_iso_mat );

if abs(Q0) > 0
    if sigma_a > 0
        add = (1-exp(-sigma_a*tau))/sigma_a * (Q0/2);
    else
        add = tau * (Q0/2);
    end
    psi = psi + add;
end
end
