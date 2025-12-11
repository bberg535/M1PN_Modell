function r = reaction(un, sigma_a, sigma_s)
    % Reaktionsterm sigma * u
    r0 = sigma_a .* un(1,:);
    r1 = (sigma_a + sigma_s) .* un(2,:);

    r = [r0;r1];
end