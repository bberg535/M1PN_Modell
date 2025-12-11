function holo_mcl_m1_pn_demo
% HOLO-MCL Kopplung von M1 und PN in EINEM File (Skeleton)
% --------------------------------------------------------
% - M1: 0. und 1. Moment U = [U0; U1]
% - PN: K Freiheitsgrade uPN(:,j) pro Zelle
% - LO-Flux: Lax-Friedrichs mit M1-Closure (Eddington)
% - HO-Flux: zentraler Flux mit PN-Closure (Eddington aus PN)
% - F_MCL = calc_g_star(HO_Flux, LO_Flux, ...)
% - Zeitdiskretisierung: SSP-RK2 (Heun)
% - Danach: Least-Squares-Fitting uPN^new auf uPN^new,HO
%           unter den Constraints M*uPN = [U0;U1]

%% Diskretisierungsparameter
L  = 2.0;
dz = 0.002;
Nz = round(L/dz);
dt = 0.001;
Tfinal = 1.0;
Nt = round(Tfinal/dt);

z = (dz/2):dz:(L - dz/2);

%% Materialparameter (einfaches Beispiel)
sigma_a = 1.0 * ones(1, Nz);   % Absorption
sigma_s = 0.0 * ones(1, Nz);   % Streuung

%% Anfangsdaten M1: Dichte-Peak in der ersten Zelle
U = zeros(2, Nz);      % [U0; U1]
U(1,1) = 1/dz;         % Masse ~1 auf erster Zelle, Rest 0
U(2,:) = 0;            % Anfangsfluss = 0

%% PN-Parameter und Initialisierung
% K = Anzahl PN-DOFs pro Zelle
K   = 4;                       % z.B. P3 -> 4 DOFs (hier nur Skeleton)
uPN = zeros(K, Nz);

% Annahme im Skeleton: erste beiden DOFs tragen U0,U1, Rest = 0
uPN(1,:) = U(1,:);
uPN(2,:) = U(2,:);
uPN(3:K,:) = 0;

% Momentmatrix M : [U0;U1] = M * uPN
% ----------------------------------------------------------
% TODO: Für dein echtes PPN: M an deine Basis anpassen!
M = [eye(2), zeros(2, K-2)];   % hier: DOF1=U0, DOF2=U1, Rest egal
MMt_inv = inv(M*M.');

%% Zeitschleife
for n = 1:Nt
    t = n*dt;

    % 1) Heun-Schritt für M1, PN nur "vorläufig" (uPN_HO)
    [U, uPN_HO] = heun_step_m1_pn(U, uPN, dt, dz, sigma_a, sigma_s);

    % 2) Least-Squares-Fitting: uPN_HO -> uPN mit Moment-Constraints U0,U1
    U_target = [U(1,:); U(2,:)];   % 2×Nz Zielmomente aus M1
    uPN = ls_fit_PN(uPN_HO, U_target, M, MMt_inv);

    % Diagnose: Gesamtmasse
    mass = sum(U(1,:)) * dz;
    fprintf('t = %.3f, mass = %.10e\n', t, mass);
end

%% Plot
U2 = calc_U2_M1(U);          % M1-Closure für 2. Moment
figure;
plot(z, U(1,:), 'k-', ...
     z, U(2,:), 'b--', ...
     z, U2,      'r:','LineWidth',1.5);
legend('U_0','U_1','U_2 (M1-Closure)','Location','best');
xlabel('z'); ylabel('Momente');
title('HOLO-MCL M1–PN Demo (Skeleton)');

end % --- Ende Hauptfunktion ---


%% ========= SSP Heun für M1 + PN-Skeleton =====================
function [U_new, uPN_HO] = heun_step_m1_pn(U, uPN, dt, dz, sigma_a, sigma_s)

Nz   = size(U,2);
im1c = [Nz, 1:Nz-1];

% ---------- Stufe 1 ----------
F1 = holo_flux_M1(U, uPN, dz, dt);          % HOLO-MCL-Fluss
q  = source_term(Nz);                      % Quelle
R1 = reaction_term(U, sigma_a, sigma_s);   % Reaktion

rhs1 = -(F1 - F1(:,im1c))/dz + q + R1;
U1   = U + dt*rhs1;

% PN: Platzhalter-Update für HO-Schritt
uPN1 = pn_step_stub(uPN, dt);              % TODO: hier dein PPN-SIU einbauen

% ---------- Stufe 2 ----------
F2 = holo_flux_M1(U1, uPN1, dz, dt);
R2 = reaction_term(U1, sigma_a, sigma_s);
rhs2 = -(F2 - F2(:,im1c))/dz + q + R2;
U2   = U1 + dt*rhs2;

% Mittelung
U_new = 0.5*(U + U2);
uPN_HO = uPN1;   % "HO"-PN-Zustand für das LS-Fitting
end


%% ========= HOLO-MCL-Fluss für M1 =============================
function F = holo_flux_M1(U, uPN, dz, dt)
% U   : 2×Nz (U0,U1)
% uPN : K×Nz
% F   : 2×Nz, begrenzter Fluss an Interfaces

Nz = size(U,2);

% --- LO-Closure (M1, Paul) ---
U2_LO = calc_U2_M1(U);                      % 1×Nz
F_LO  = lo_flux_M1(U, U2_LO);               % 2×Nz

% --- HO-Closure aus PN ---
U2_HO = calc_U2_from_PN(uPN);               % 1×Nz (TODO: PPN-Formel!)
F_HO  = ho_flux_M1(U, U2_HO);               % 2×Nz

% --- MCL-Korrektur (FCT / G^*) ---
flux_fun    = @(Ustate) m1_flux(Ustate, calc_U2_M1(Ustate));
wavespeed   = @(Ustate) ones(1,size(Ustate,2));   % |λ|max ≈ 1
F = calc_g_star(F_HO, F_LO, flux_fun, wavespeed, Nz, U);
end


%% ========= LO-Flux (Lax-Friedrichs) ==========================
function F_LO = lo_flux_M1(U, U2)
% U  : 2×Nz, U2: 1×Nz
[nVar, Nz] = size(U);
ip0c = 1:Nz;
ip1c = [2:Nz 1];

ul = U(:,ip0c);             % U_i
ur = U(:,ip1c);             % U_{i+1}
U2L = U2(ip0c);
U2R = U2(ip1c);

FL = m1_flux(ul, U2L);
FR = m1_flux(ur, U2R);

lambda    = 1;                         % für M1 typ. ≤1 (c=1)
lambdaMat = lambda * ones(nVar, Nz);

F_LO = 0.5*(FL + FR) - 0.5*lambdaMat .* (ur - ul);
end


%% ========= HO-Flux (zentraler Flux) ==========================
function F_HO = ho_flux_M1(U, U2)
% zentraler (nicht diffusativer) Flux
Nz = size(U,2);
ip0c = 1:Nz;
ip1c = [2:Nz 1];

ul = U(:,ip0c);
ur = U(:,ip1c);
U2L = U2(ip0c);
U2R = U2(ip1c);

FL = m1_flux(ul, U2L);
FR = m1_flux(ur, U2R);

F_HO = 0.5*(FL + FR);
end


%% ========= Physikalischer Fluss des M1-Systems ===============
function F = m1_flux(U, U2)
% U:  2×Nz  (U0,U1)
% U2: 1×Nz  (2. Moment)
F = [U(2,:); U2];
end


%% ========= M1-Closure: 2. Moment =============================
function U2 = calc_U2_M1(U)
E = U(1,:);
F1 = U(2,:);
eps_rho = 1e-14;

f  = abs(F1) ./ (E + eps_rho);              % Flux-Faktor
chi = (3 + 4*f.^2) ./ (5 + 2*sqrt(4 - 3*f.^2));  % Eddington-Faktor
U2 = chi .* E;
end


%% ========= 2. Moment aus PN (HO-Closure) =====================
function U2 = calc_U2_from_PN(uPN)
% Skeleton:
% --------------------------------------------
% TODO: HIER deine Positive-PN-Formel einsetzen:
%       aus uPN die Winkelverteilung rekonstruieren
%       und U2^HO = ∫ μ^2 ψ_N(uPN,μ) dμ berechnen.
%
% Im Skeleton nehmen wir einfach an:
%  - DOF1 ≈ U0, DOF2 ≈ U1
%  - und benutzen denselben M1-Closure.
U0 = uPN(1,:);
U1 = uPN(2,:);
U  = [U0; U1];
U2 = calc_U2_M1(U);
end


%% ========= Reaktions-Term (Absorption/Streuung) ==============
function R = reaction_term(U, sigma_a, sigma_s)
E = U(1,:);
F1 = U(2,:);
R0 = - sigma_a .* E;
R1 = - (sigma_a + sigma_s) .* F1;
R  = [R0; R1];
end


%% ========= Quellterm (einfaches Beispiel) ====================
function q = source_term(Nz)
q0 = zeros(1, Nz);
q1 = zeros(1, Nz);

mask = zeros(1, Nz);
mask(1) = 1;          % Quelle nur in der ersten Zelle

S0  = 1.0;            % Stärke
q0  = S0 * mask;

f_src = 0.0;          % Anisotropie der Quelle
q1  = f_src .* q0;

q   = [q0; q1];
end


%% ========= PN-Update-Platzhalter =============================
function uPN_new = pn_step_stub(uPN, dt)
% TODO: Hier deinen echten PPN-SIU-Step einbauen.
%       Für das Skeleton lassen wir PN einfach konstant.
uPN_new = uPN;
end


%% ========= Least-Squares-Fitting PN ==========================
function uPN = ls_fit_PN(uPN_HO, U_target, M, MMt_inv)
% uPN_HO  : K×Nz  (PN nach HO-Schritt)
% U_target: 2×Nz  (Zielmomente [U0;U1] aus M1)
% M       : 2×K   Momentmatrix
% MMt_inv : (M*M')^{-1}, 2×2

[K, Nz] = size(uPN_HO);
uPN = zeros(K, Nz);

for j = 1:Nz
    vHO  = uPN_HO(:,j);
    Utar = U_target(:,j);

    correction = M' * (MMt_inv * (Utar - M*vHO));
    uPN(:,j)   = vHO + correction;
end
end


%% ========= MCL-Flux-Limiter (G*) =============================
function Gstar = calc_g_star(HO_Flux, LO_Flux, flux, fluxj, Nz, u)
% HO_Flux, LO_Flux : nVar×Nz, Interfaceflüsse
% flux             : @(u) -> nVar×Nz, physikalischer Fluss
% fluxj            : @(u) -> 1×Nz, max. Wellengeschwindigkeit
% u                : nVar×Nz, Zellmittelwerte
[nVar, Nz_check] = size(u);
if Nz_check ~= Nz
    error('calc_g_star: size(u,2) ~= Nz');
end

ip0c = 1:Nz;
ip1c = [2:Nz 1];
im1c = [Nz 1:Nz-1];

ulc = u(:,ip0c);              % u_i
urc = u(:,ip1c);              % u_{i+1}

lamL = abs(fluxj(ulc));       % 1×Nz
lamR = abs(fluxj(urc));       % 1×Nz
lambda = max(lamL, lamR);     % 1×Nz
lambdaMat = repmat(lambda, nVar, 1);

% antidiffusiver Anteil
fAe = HO_Flux - LO_Flux;      % nVar×Nz

% lokale Min/Max-Schranken
umax = max(u(:,im1c), max(u(:,ip0c), u(:,ip1c)));
umin = min(u(:,im1c), min(u(:,ip0c), u(:,ip1c)));

% Bar-State (Rusanov)
F_L = flux(ulc);
F_R = flux(urc);
wbar = 0.5*(urc + ulc) - 0.5*(F_R - F_L) ./ lambdaMat;

% begrenzter antidiffusiver Flux (MCL-Formel)
fAe = ...
    min( max(0, fAe), ...
         min(lambdaMat.*umax(:,ip0c) - wbar, ...
             wbar - lambdaMat.*umin(:,ip1c)) ) ...
  + max( min(0, fAe), ...
         max(lambdaMat.*umin(:,ip0c) - wbar, ...
             wbar - lambdaMat.*umax(:,ip1c)) );

% finaler Fluss
Gstar = LO_Flux + fAe;
end
