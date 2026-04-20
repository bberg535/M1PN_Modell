N = [1,3];

% Quadraturpunkte
[omega, w] = gausslegendre(10);

% Referenzlösung
f = @(x) (0.1./(exp(4.*x)));       % Schreibfehler im Paper (/ statt *)
f_omega = f(omega);

figure;
for i = 1:2
    % Auswertung der Legendrepolynome an den Quadraturpunkten
    p = legpoly_eval(omega,N(i));

    % Quadratur von u; s. Abschnitt unter (2)
    u = p * (w .* f_omega);

    % Optimierungsproblem ergibt PN Rekonstruktion
    PN = quadprog(diag(w),zeros(size(omega)),[],[],p*diag(w),u,[],[]);

    % Hinzugefügter Lower-Bound
    PPN = quadprog(diag(w),zeros(size(omega)),[],[],p*diag(w),u,zeros(size(omega)),[]);

    subplot(1,2,i)
    hold on
    xlabel('\mu'); title('P_')
    scatter(omega, f_omega, 36, 'o', 'filled');
    scatter(omega, PN,      60, '*');
    scatter(omega, PPN,     36, 'x', 'LineWidth',1.2);
    title(sprintf('(%c) P%d vs PP%d', 'a'+i-1, N(i), N(i))); 
    hold off
end


%% Hilfsfunktionen

function [x,w] = gausslegendre(n)
    i = (1:n-1).';
    beta = i ./ sqrt(4*i.^2 - 1);
    J = diag(zeros(n,1)) + diag(beta,1) + diag(beta,-1);
    [V,D] = eig(J);
    [x,idx] = sort(diag(D));
    V = V(:,idx);
    w = 2*(V(1,:)'.^2);
end

function P = legpoly_eval(mu, N)
    mu = mu(:).'; M = numel(mu);
    P = zeros(N+1, M);
    P(1,:) = 1;
    if N>=1, P(2,:) = mu; end
    for l = 1:N-1
        P(l+2,:) = ((2*l+1)*mu.*P(l+1,:) - l*P(l,:)) / (l+1);
    end
end
