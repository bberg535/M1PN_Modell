function u = heuns_method(un, dt, dx, f, fj, Nz, eddington, method)
    ip0c = (1:Nz)'; ip1c = [2:Nz 1]';

    ulc = un(:,ip0c); urc = un(:,ip1c);
    lambda = max(abs(fj(ulc)),abs(fj(urc)));

    % Stufe 1
    eps_rho = 1e-14;
    f0   = abs(un(2,:))./(un(1,:) + eps_rho);
    psi2_0 = eddington(f0).*un(1,:);
    f1  = method(un, dt, dx, f, lambda, psi2_0, Nz);
    u1  = un + dt .* f1;

    % Stufe 2: psi2 aus u1!
    f1s  = abs(u1(2,:))./(u1(1,:) + eps_rho);
    psi2_1 = eddington(f1s).*u1(1,:);
    f2  = method(u1, dt, dx, f, lambda, psi2_1, Nz);
    u2 = u1 + dt .* f2;

    % Heun-Update
    %u   = un + 0.5 .* dt .* (f1 + f2);
    u = (un + u2)/2;
end