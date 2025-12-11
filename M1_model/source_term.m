function q = source_term(Nz)
    % Quelle strahlt rechts rein
    bound = zeros(1, Nz);
    bound(1) = 1;

    % Stärke der Quelle
    S0 = 10000.0;
    q0 = S0 * bound;

    % Richtungsparameter f in [-1, 1], sodass abs(q1) <= q0
    f_src = 0.0;
    q1 = f_src .* q0;

    q = [q0; q1];
end