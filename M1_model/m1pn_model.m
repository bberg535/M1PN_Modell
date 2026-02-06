% M1PN-Modell für die lineare Boltzmann Gleichung
% mit periodischen Randwerten und drei verschiedenen Verfahren
%
% Quelle: 
%   https://arxiv.org/abs/2509.07689v1 für ganzzahlige Gleichungen
%   https://doi.org/10.1142/13466 für alle anderen Gleichungen

clear;

%% Eingabewerte
% Verfahren bestimmen: 1 [MCL]; -1 [LLF]; sonst [LW]
ilim = 1;
cont_plot = false;

% Größen für die Diskretisierung
N = 1;         % N fest, da M1-Modell
NPN = 3;
dt = 0.001;
dz = 0.002;    % CFL: dt < dz
nec = 2/dz;
nne = nec + 1;
cfl = dt/dz;
T = 1;

z = (0 : dz : 2);
mu = z - 1;
legpols = legpoly_eval(mu,NPN);

% Flussfunktion
flux = @(u) ([u(2,:); u(3,:)]); fluxj = @(u) 1;

% Absorptions- und Streuungsterm
sigma_a = (dz .* ones(nne,1))';
sigma_s = (dz .* ones(nne,1))';
source_strength = 1.0;

% Anfangswerte
rho0 = [1/dz;zeros(nne-1,1)];
u = zeros(N+1, nne);
uPN = zeros(NPN+1, nne);
u(1,:) = rho0.'; 
uPN(1,:) = rho0.';

if cont_plot
    psi2 = calc_psi2(u);
    figure;
    plot(z,[u;psi2])
end

%% Zeitschleife
for t = 0:dt:T
    % Zeitschritt mit SSP Heuns Method
    [u, uPN] = heuns_method(u, uPN, dt, dz, flux, nec, sigma_a, sigma_s, ilim , NPN, source_strength, legpols);

    % Least Squares Fitting

    % Norm- und Positivitätsüberprüfung
    mass = sum(u(1,:)) * dz; 
    fprintf('t = %.3f, mass = %.10e\n', t, mass);
    [min_pn, idx] = min(uPN(1,:));
    fprintf('min rho_M1 = %.3e, min rho_PN = %.3e at %d\n ', ...
        min(u(1,:)), min_pn, idx);

    % Kontinuierlicher Plot
    if cont_plot
        psi2 = calc_psi2(u);
        plot(z,[u;psi2])
        %plot(z,uPN)
        drawnow;
    end
end

%% Plotten
psi2 = calc_psi2(u);
figure;
plot(z,[u;psi2])
figure;
plot(z,uPN)


