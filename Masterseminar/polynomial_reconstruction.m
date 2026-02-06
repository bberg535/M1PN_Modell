function [uplus, uminus] = polynomial_reconstruction(u, ub)
[nM,Nx] = size(u);

u_ext = [ub, u, ub];                 % (nM x (Nx+2))

dF = u_ext(:,3:end)   - u_ext(:,2:end-1);     % u_{i+1}-u_i
dB = u_ext(:,2:end-1) - u_ext(:,1:end-2);     % u_i-u_{i-1}
dC = 0.5*(u_ext(:,3:end) - u_ext(:,1:end-2)); % 1/2(u_{i+1}-u_{i-1})

s  = minmod(dF,dB,dC);              % (nM x Nx) = Δz * u'_i

uLcell = u - 0.5*s;                  % u_i(z_{i-1/2})
uRcell = u + 0.5*s;                  % u_i(z_{i+1/2})

uminus = zeros(nM, Nx+1);            % linke Zustände u^-_{i+1/2}
uplus  = zeros(nM, Nx+1);            % rechte Zustände u^+_{i+1/2}

uminus(:,1) = ub;        uplus(:,1)  = uLcell(:,1);
for i=1:Nx-1
    uminus(:,i+1) = uRcell(:,i);
    uplus(:,i+1)  = uLcell(:,i+1);
end
uminus(:,Nx+1) = uRcell(:,Nx);
uplus(:,Nx+1)  = ub;
end
