function Gstar = calc_g_star(HO_Flux, LO_Flux, f, fj, Nz, un, psi2) 
ip0c = [1:Nz]'; ip1c = [2:Nz 1]'; im1c = [Nz 1:Nz-1]'; 
ulc = un(:,ip0c); urc = un(:,ip1c); 

lambda = max(abs(fj(ulc)),abs(fj(urc))); 

fAe = HO_Flux - LO_Flux; 

% Lokale Schranken Gleichung (22) 
umax = max(un(:,im1c),max(un(:,ip0c),un(:,ip1c))); 
umin = min(un(:,im1c),min(un(:,ip0c),un(:,ip1c))); 

% Bar-State berechnen Gleichung (5) 
wbar = 0.5*(urc+ulc).*lambda-0.5*(f([urc; psi2(ip1c)])-f([ulc; psi2(ip0c)])); 

% Anteil des maximalen HO-Flusses bestimmen 
fAe = min(max(0,fAe),min(lambda.*umax(:,ip0c)-wbar,wbar-lambda.*umin(:,ip1c))) ...
    +max(min(0,fAe),max(lambda.*umin(:,ip0c)-wbar,wbar-lambda.*umax(:,ip1c)));
% Vom ursprünglichen Fluss abziehen 
Gstar = LO_Flux - fAe; 

end