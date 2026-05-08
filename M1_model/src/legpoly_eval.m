function P = legpoly_eval(mu,N)
% Berechnung der Legendre-Polynome durch die Bonnet-Formel
    mu = mu(:).';
    M = numel(mu);
    P = zeros(N+1,M);
    P(1,:) = 1;
    if N>=1, P(2,:) = mu; end
    for l = 1:N-1
        P(l+2,:) = ((2*l+1)*mu.*P(l+1,:) - l*P(l,:))/(l+1);
    end
end