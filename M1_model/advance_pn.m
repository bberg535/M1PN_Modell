function u = advance_pn(un, N, method, dt, dz)
Nz = size(un,2); cfl = dt / dz;
ip0c = [1:Nz]'; ip1c = [2:Nz 1]'; im1c = [Nz 1:Nz-1]';



% Transformations- und Basisvektoren/Matrizen
[xi,w] = gausslegendre(N+1);
P = legpoly_eval(xi, N);
G = diag(2./(2.*(0:N)+1));
R = P*diag(w);
L = (G \ P).';

flux = @(v) (xi.*v); fluxj = @(v) xi;

v = L * un;

vlc=v(:,ip0c); vrc=v(:,ip1c); ve=0.5*(vlc+vrc);

lambdaMCL=max(abs(fluxj(vlc)),abs(fluxj(vrc)));

sigma = 1;
gamma = 1/(1+sigma*dt);

% Central difference flux
fCD=0.5*(flux(vrc)+flux(vlc));

% Local Lax-Friedrichs flux 
fLF=fCD-0.5*lambdaMCL.*(vrc-vlc); % Gleichung (4)

% Lax-Wendroff flux   
fLW=fCD-0.5*cfl*(lambdaMCL.^2).*(vrc-vlc);

if method == 1

    % Antidiffusive flux
    fAe=fLF-fLW; % Gleichung (11)

    % Local bounds
    umax=max(v(:,im1c),max(v(:,ip0c),v(:,ip1c)));
    umin=min(v(:,im1c),min(v(:,ip0c),v(:,ip1c)));

    % Scaled bar states
    wbar=0.5*(vrc+vlc).*lambdaMCL-0.5*(flux(vrc)-flux(vlc));

    % Flux limiting
    fAe=min(max(0,fAe),min(lambdaMCL.*umax(:,ip0c)-wbar,wbar-lambdaMCL.*umin(:,ip1c))) ...
        +max(min(0,fAe),max(lambdaMCL.*umin(:,ip0c)-wbar,wbar-lambdaMCL.*umax(:,ip1c)));

    % Flux correction
    fMCL=fLF-fAe; % Gleichung (18)

elseif method == -1

    fMCL=fLF;      

else

    fMCL=fLW;

end

rho = w.'*(v(:,ip0c)-cfl*(fMCL(:,ip0c)-fMCL(:,im1c)));

v(:,ip0c)=gamma*(v(:,ip0c)-cfl*(fMCL(:,ip0c)-fMCL(:,im1c))) + (1-gamma)*rho/2;

u = R * v;
end