function u = MCL(un, dt, dx, f, lambda, psi2, Nz)
    ip0c = [1:Nz]'; ip1c = [2:Nz 1]'; im1c = [Nz 1:Nz-1]';
    ulc = un(:,ip0c); urc = un(:,ip1c);
    fCD = 0.5*(f([urc; psi2(ip1c)])+f([ulc; psi2(ip0c)]));
    fLF = fCD - 0.5 * lambda .* (urc - ulc);
    fLW = fCD - 0.5 * dt/dx * lambda.^2 .* (urc - ulc);
    fAe = fLF - fLW;
    umax = max(un(:,im1c),max(un(:,ip0c),un(:,ip1c)));
    umin = min(un(:,im1c),min(un(:,ip0c),un(:,ip1c)));
    wbar = 0.5*(urc+ulc).*lambda-0.5*(f([urc; psi2(ip1c)])-f([ulc; psi2(ip0c)]));
    fAe = min(max(0,fAe),min(lambda.*umax(:,ip0c)-wbar,wbar-lambda.*umin(:,ip1c))) ...
            +max(min(0,fAe),max(lambda.*umin(:,ip0c)-wbar,wbar-lambda.*umax(:,ip1c)));
    u = fLF - fAe;
end
