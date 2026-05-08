function u = advance_pn(un, N, method, dt, dz)
% Berechnung der Rechten Seite des PN-Modells aufgebaut nach 
%
% https://doi.org/doi/10.1137/15M1052871
%
% Mit FV-Verfahren LLF, LW oder MCL

Nz = size(un,2) - 1; cfl = dt / dz;
ip0c = [1:Nz]'; ip1c = [2:Nz 1]'; im1c = [Nz 1:Nz-1]';
rhs_v = zeros(N+1, Nz + 1);

% Transformations- und Basisvektoren/Matrizen
[xi,w] = gausslegendre(N+1);
P = legpoly_eval(xi, N);
G = diag(2./(2.*(0:N)+1));
R = P*diag(w);
L = (G \ P).';

% Diagonalisiere das System per Basistransformation
v = L * un;
vlc = v(:,ip0c); vrc = v(:,ip1c);

% Flussfunktion und -geschwindigkeit
flux = @(v) (xi.*v); fluxj = @(v) xi;
lambdaMCL=max(abs(fluxj(vlc)),abs(fluxj(vrc)));

% Central difference flux
fCD=0.5*(flux(vrc)+flux(vlc));

% Local Lax-Friedrichs flux 
fLF=fCD-0.5*lambdaMCL.*(vrc-vlc);

% Lax-Wendroff flux   
fLW=fCD-0.5*cfl*(lambdaMCL.^2).*(vrc-vlc);

if method == 1
    fAe=fLF-fLW;
    umax=max(v(:,im1c),max(v(:,ip0c),v(:,ip1c)));
    umin=min(v(:,im1c),min(v(:,ip0c),v(:,ip1c)));
    wbar=0.5*(vrc+vlc).*lambdaMCL-0.5*(flux(vrc)-flux(vlc));
    fAe=min(max(0,fAe),min(lambdaMCL.*umax(:,ip0c)-wbar,wbar-lambdaMCL.*umin(:,ip1c))) ...
        +max(min(0,fAe),max(lambdaMCL.*umin(:,ip0c)-wbar,wbar-lambdaMCL.*umax(:,ip1c)));
    fMCL=fLF-fAe;
elseif method == -1
    fMCL=fLF;      
else
    fMCL=fLW;
end

rhs_v(:,ip0c) = -(fMCL(:,ip0c) - fMCL(:,im1c)) / dz;
u = R * rhs_v;  
end
