function [x, w] = gauss_legendre(n, a, b)
% Golub-Welsch Gauss-Legendre quadrature on [a,b]
% Returns nodes x and weights w.
i = (1:n-1)';
beta = i ./ sqrt(4*i.^2 - 1);
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x = diag(D);
[w, idx] = sort(x);
V = V(:,idx);
x = w;
w = 2*(V(1,:).^2)';

% map from [-1,1] to [a,b]
x = (b-a)/2 * x + (a+b)/2;
w = (b-a)/2 * w;
end
