function P = legendreP_matrix(N, mu)
% Returns P_l(mu) for l=0..N as matrix (N+1 x numel(mu))
mu = mu(:).';
P = zeros(N+1, numel(mu));
P(1,:) = 1;
if N==0, return; end
P(2,:) = mu;
for l=1:N-1
    % recurrence: (l+1)P_{l+1} = (2l+1)mu P_l - l P_{l-1}
    P(l+2,:) = ((2*l+1).*mu.*P(l+1,:) - l*P(l,:)) / (l+1);
end
end