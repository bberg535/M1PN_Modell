% Messung
tAll = tic;

N = [1,3,5,7];
pn_img = cell(numel(N),3);
ppn_img = cell(numel(N),3);
for j = 1 : numel(N)
    % Größen für die Diskretisierung
    M = 2*(N(j)+1);
    dt = 0.01;
    dz = 0.04;    % CFL: dt < dz
    L = 24;
    Nz = round(L/dz);
    index = floor(Nz/2) + 1;
    z  = (-L/2 + dz/2) : dz : (L/2 - dz/2);
    
    % Quadraturpunkt und -gewichte berechnen
    [omega,w] = gausslegendre(M);
    P = legpoly_eval(omega, N(j));
    pos = omega > 0;

    % Faktoren 
    sigma = 1;
    Q = diag([0; ones(N(j),1)]);
    G = diag(2./(2.*(0:N(j))+1));
    gamma = eye(N(j)+1)/(eye(N(j)+1) + dt*sigma*Q);

    % Anfangswerte
    u_pn = zeros(N(j)+1, Nz);
    u_pn(1,index) = 1/dz;
    u_ppn = u_pn;

    % Zeitschleife
    for t = 0:dt:5
        % Normale PN-Berechnung
        c_pn = G \ u_pn;
        F_PN = P.' * c_pn;

        c_ppn = G \ u_ppn;
        F_PPN = P.' * c_ppn;

        % Abfrage auf negative Werte in F
        neg = any(F_PPN < 0, 1);
        F = F_PPN;
        
        % Erhalte Ansatz aus Optimierungsproblem 
        if any(neg) 
            for i = find(neg) 
               F(:,i)  = quadprog(diag(w), zeros(size(omega)), [], [], P*diag(w), u_ppn(:,i), zeros(size(omega)), []); 
            end 
        end

        % Indexshift für periodische Randwerte
        sl = [Nz, 1:Nz-1];
        sr = [2:Nz, 1];

        % SIU (33)
        Spos_pn = P(:,pos) * (w(pos).*omega(pos).*(F_PN(pos,:) - F_PN(pos,sl))/dz);
        Smin_pn = P(:,~pos) * (w(~pos).*omega(~pos).*(F_PN(~pos,sr) - F_PN(~pos,:))/dz);
        u_pn = gamma * (u_pn - dt* (Spos_pn + Smin_pn));

        Spos_ppn = P(:,pos) * (w(pos).*omega(pos).*(F(pos,:) - F(pos,sl))/dz);
        Smin_ppn = P(:,~pos) * (w(~pos).*omega(~pos).*(F(~pos,sr) - F(~pos,:))/dz);
        u_ppn = gamma * (u_ppn - dt* (Spos_ppn + Smin_ppn));

        if t == 1
            pn_img{j,1} = u_pn(1,:);
            ppn_img{j,1} = u_ppn(1,:);
        end
        if t == 2
            pn_img{j,2} = u_pn(1,:);
            ppn_img{j,2} = u_ppn(1,:);
        end
        if t == 5
            pn_img{j,3} = u_pn(1,:);
            ppn_img{j,3} = u_ppn(1,:);
        end
    end
end

tEnd = toc(tAll);

%% Plots

window = {(z>=-2 & z<= 2),(z>=-3 & z<= 3),(z>=-6 & z<= 6)};
titles = ["t = 1.0","t = 2.0","t = 5.0"];
for h = 1:3
    figure;
    for j = 1:numel(N)
        subplot(2,2,j); hold on
        plot(z(window{h}), pn_img{j,h}(window{h}), 'o', 'MarkerSize', 3);
        plot(z(window{h}), ppn_img{j,h}(window{h}), '+', 'MarkerSize', 3);
        title(sprintf('P%d–S%d', N(j), 2*(N(j)+1)));
        hold off
    end
    sgtitle(sprintf('Plane source, %s', titles(h)));
end


%% Hilfsfunktionen

function [x,w] = gausslegendre(n)
    i = (1:n-1)';
    beta = i./sqrt(4*i.^2 - 1);
    J = diag(zeros(n,1)) + diag(beta,1) + diag(beta,-1);
    [V,D] = eig(J);
    [x,idx] = sort(diag(D));
    V = V(:,idx);
    w = 2*(V(1,:)'.^2);
end

function P = legpoly_eval(mu,N)
    mu = mu(:).';
    M = numel(mu);
    P = zeros(N+1,M);
    P(1,:) = 1;
    if N>=1, P(2,:) = mu; end
    for l = 1:N-1
        P(l+2,:) = ((2*l+1)*mu.*P(l+1,:) - l*P(l,:))/(l+1);
    end
end