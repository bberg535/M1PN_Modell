% Low-order scheme and MCL for the 1D advection equation 
% 
%      u_t+(f(u))_x=0 with f(u)=u
% 
% in (xl,xr) with periodic boundary conditions

% Reference: Section 3, equation for h in
%
% https://doi.org/10.1016/j.matcom.2025.04.031

clear

% Boundary points

xl=0; xr=1;

% Limiting on (ilim=1), off (ilim=0), low-order (ilim=-1)

ilim=1;

% Discretization parameters

cfl=0.5; % constant factor dt/dx 

nec=200; % number of mesh elements

% Tolerance
eps=1e-15;

% Mesh data
dxc=(xr-xl)/nec;           % mesh size
nnc=nec+1;                 % number of nodes
xc=linspace(xl,xr,nnc);    % coordinates of nodes
dtc=cfl*dxc;               % time step

% Final time T, initial data, flux f(u), Jacobian f'(u)

T=0.5; uc=zeros(nnc,1);

uc(0.1 <= xc & xc <= 0.3)=1;
    
flux = @(u) u; fluxj = @(u) 1;

% Index arrays (i+0, i+1, i-1)
ip0c=[1:nec]'; ip1c=[2:nec 1]'; im1c=[nec 1:nec-1]';

% Time stepping loop

t=0; 

while t+dtc <= T+eps

    t = t+dtc;

    % Left and right Riemann states
    ulc=uc(ip0c); urc=uc(ip1c); ue=0.5*(ulc+urc);

    % Rusanov maximum speed 
    lambda=max(abs(fluxj(ulc)),abs(fluxj(urc))); % Gleichung (6)
    
    % Central difference flux
    fCD=0.5*(flux(urc)+flux(ulc));

    % Local Lax-Friedrichs flux 
    fLF=fCD-0.5*lambda.*(urc-ulc); % Gleichung (4)

    % Lax-Wendroff flux   
    fLW=fCD-0.5*cfl*(lambda.^2).*(urc-ulc); % LW korrespondiert zu feiner Gitter Approximation

    if ilim == 1

      % Antidiffusive flux
      fAe=fLF-fLW; % Gleichung (11)

      % Local bounds
      umax=max(uc(im1c),max(uc(ip0c),uc(ip1c))); % Gleichung (24a) mit Gleichung (22)
      umin=min(uc(im1c),min(uc(ip0c),uc(ip1c))); % Gleichung (24a) mit Gleichung (22) 

      % Scaled bar states
      % Im Buch (3.72) aber mit lambda = 2*d_ij multipliziert
      wbar=0.5*(urc+ulc).*lambda-0.5*(flux(urc)-flux(ulc)); 

      % Flux limiting
      %                  % f_ij max (Unter Gleichung (3.77))
      fAe=min(max(0,fAe),min(lambda.*umax(ip0c)-wbar,wbar-lambda.*umin(ip1c))) ...
         +max(min(0,fAe),max(lambda.*umin(ip0c)-wbar,wbar-lambda.*umax(ip1c))); % Gleichung (27)

      % Flux correction
      fMCL=fLF-fAe; % Gleichung (18)

    elseif ilim == -1

      fMCL=fLF;      

    else

      fMCL=fLW;

    end
    
    % Explicit update
    uc(ip0c)=uc(ip0c)-cfl*(fMCL(ip0c)-fMCL(im1c)); % Gleichung (20)

    % Periodic boundary conditions
    uc(nnc)=uc(1); 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); plot(xc,uc,'LineWidth',2);
title('Numerical solution')
