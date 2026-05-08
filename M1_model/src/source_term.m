function q = source_term(Nz, source_strength)
    % Quelle strahlt rechts rein
    bound = zeros(1, Nz);
    bound(1) = 1;

    % Stärke der Quelle
    q0 = source_strength * bound;

    % Richtungsparameter f in [-1, 1], sodass abs(q1) <= q0
    f_src = 0.0;
    q1 = f_src .* q0;

    % Endergebnis
    q = [q0; q1];
end