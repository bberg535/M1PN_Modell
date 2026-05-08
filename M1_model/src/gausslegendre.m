function [x,w] = gausslegendre(n)
% Gauß-Legendre Quadratur nach dem Golub-Welsch Algorithmus
    i = (1:n-1)';
    beta = i./sqrt(4*i.^2 - 1);
    J = diag(zeros(n,1)) + diag(beta,1) + diag(beta,-1);
    [V,D] = eig(J);
    [x,idx] = sort(diag(D));
    V = V(:,idx);
    w = 2*(V(1,:)'.^2);
end