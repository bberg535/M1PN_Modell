% M1-Modell für die lineare Boltzmann Gleichung
% mit periodischen Randwerten und drei verschiedenen Verfahren
%
% Quelle: 
%   https://arxiv.org/abs/2509.07689v1 für ganzzahlige Gleichungen
%   https://doi.org/10.1142/13466 für alle anderen Gleichungen

clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(script_dir)));
paths = setup_project_paths(); %#ok<NASGU>

%% Eingabewerte
% Verfahren bestimmen: 1 [MCL]; -1 [LLF]; sonst [LW]
ilim = 1;
method = select_flux_method(ilim);

% Größen für die Diskretisierung
N = 1;         % N fest, da M1-Modell
dt = 0.001;
dz = 0.002;    % CFL: dt < dz
Nz = 2/dz;
cfl = dt/dz;
T = 1.5;

% Flussfunktion
flux = @(u) ([u(2,:); u(3,:)]); fluxj = @(u) 1;

% Absorptions- und Streuungsterm
sigma_a = 0 * (dz .* ones(Nz,1))';
sigma_s = 0 * (dz .* ones(Nz,1))';
source_strength = 0.0;

% Anfangswerte
rho0 = [1/dz;zeros(Nz-1,1)];
u = zeros(N+1, Nz);
u(1,:) = rho0.'; 



%% Zeitschleife
num_steps = round(T / dt);
if abs(num_steps * dt - T) > 1e-12 * max(1, T)
    error('m1_model:timeGridMismatch', ...
        'T=%.16g is not an integer multiple of dt=%.16g.', T, dt);
end

for step = 1:num_steps
    t = step * dt;
    % Zeitschritt mit SSP Heuns Method
    u = advance_m1(u, dt, dz, flux, fluxj, Nz, sigma_a, sigma_s, method, source_strength);
    mass = sum(u(1,:)) * dz; 
    fprintf('t = %.3f, mass = %.10e\n', t, mass);
end

%% Plotten
psi2 = calc_psi2(u);
z = (dz/2 : dz : 2 - dz/2);

figure;
plot(z,[u;psi2])
