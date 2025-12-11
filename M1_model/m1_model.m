% M1-Modell für die lineare Boltzmann Gleichung
% mit periodischen Randwerten und drei verschiedenen Verfahren
%
% Quelle: 
%   https://arxiv.org/abs/2509.07689v1 für ganzzahlige Gleichungen
%   https://doi.org/10.1142/13466 für alle anderen Gleichungen

clear;

%% Eingabewerte
% Verfahren bestimmen: 1 [MCL]; -1 [LLF]; sonst [LW]
ilim = 1;

% Größen für die Diskretisierung
N = 1;         % N fest, da M1-Modell
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
u(1,:) = rho0.'; 



%% Zeitschleife
for t = 0:dt:T
    switch ilim
        case 1
            method = @MCL;
        case -1
            method = @LLF;
        otherwise
            method = @LW;
    end  
    
    % Zeitschritt mit SSP Heuns Method
    u = advance_m1(u, dt, dz, flux, fluxj, Nz, sigma_a, sigma_s, method);
    mass = sum(u(1,:)) * dz; 
    fprintf('t = %.3f, mass = %.10e\n', t, mass);
end

%% Plotten
psi2 = calc_psi2(u);
z = (dz/2 : dz : 2 - dz/2);

figure;
plot(z,[u;psi2])
