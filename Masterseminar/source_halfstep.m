function u = source_halfstep(u, tau, sigma_s, sigma_a, Q0)
% Exact update for
% d/dt psi = -sigma_t psi + sigma_s * (phi/2) + Q0/2 (isotropic source)
% In moment form: u0 reacts to Q0/absorption, higher moments damp with exp(-sigma_t t)
%
% For the plane-source test typically: sigma_a=0, Q0=0 -> u0 constant, higher moments damp exp(-sigma_s*t).

sigma_t = sigma_s + sigma_a;

% map u to isotropic state 
u_iso = u;
u_iso(2:end,:) = 0;

% Gleichung (4.12)
ea = exp(-sigma_a*tau);
es = exp(-sigma_s*tau);

u = ea*( es*u + (1-es)*u_iso );

if abs(Q0) > 0
    % (b) ist für Legendre-Basis
    % (2,0,0,0,...), daher wird nur auf u0 addiert
    if sigma_a > 0
        u(1,:) = u(1,:) + (1-exp(-sigma_a*tau))/sigma_a * Q0;
    else
        % Grenzwertbildung für sigma a gegen 0, tau = dt/2
        u(1,:) = u(1,:) + tau * Q0;
    end
end
end