% M1-Modell für die lineare Boltzmann Gleichung
% mit periodischen Randwerten und drei verschiedenen Verfahren
%
% Quelle: 
%   https://arxiv.org/abs/2509.07689v1 für ganzzahlige Gleichungen
%   https://doi.org/10.1142/13466 für alle anderen Gleichungen

clear;

%% Eingabewerte
% Verfahren bestimmen: 1 [MCL]; -1 [LLF]; sonst [LW]
ilim = -1;

% Größen für die Diskretisierung
N = 1;         % N fest, da M1-Modell
NPN = 5;
dt = 0.001;
dz = 0.002;    % CFL: dt < dz
Nz = 2/dz;
cfl = dt/dz;
T = 1;

% Flussfunktion
flux = @(u) ([u(2,:); u(3,:)]); fluxj = @(u) 1;

% Absorptions- und Streuungsterm
sigma_a = (dz .* ones(Nz,1))';
sigma_s = (dz .* ones(Nz,1))';

% Anfangswerte
rho0 = [1/dz;zeros(Nz-2,1); 1/dz];
u = zeros(N+1, Nz);
uPN = zeros(NPN+1, Nz);
u(1,:) = rho0.'; 
uPN(1,:) = rho0.';

%% Zeitschleife
for t = 0:dt:T
    % Zeitschritt mit SSP Heuns Method
    [u, uPN] = heuns_method(u, uPN, dt, dz, flux, fluxj, Nz, sigma_a, sigma_s, ilim , NPN);

    % U_new = [u(1,:); u(2,:)];

    psi2HO = (uPN(1,:) + 2 .* uPN(3,:)) / 3;

    uHO = [uPN(1,:); uPN(2,:)];

    uPN = ls_fit_PN_legendre(uPN, u);

    fprintf('min rho_M1 = %.3e, min rho_PN = %.3e\n', ...
        min(u(1,:)), min(uPN(1,:)));
    mass = sum(u(1,:)) * dz; 
    fprintf('t = %.3f, mass = %.10e\n', t, mass);
end

%% Plotten
psi2 = calc_psi2(u);
z = (dz/2 : dz : 2 - dz/2);

figure;
plot(z,[u;psi2])
figure;
plot(z,uPN)
