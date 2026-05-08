function Gstar = real_HO_flux(uPN, u, dt, dx, ilim, flux, nec, ~, ~, ~, ~)
% Berechne HO- und LO-Interface-Fluesse und limitiere den Korrekturfluss.

    nne = nec + 1;
    if size(u, 2) ~= nne || size(uPN, 2) ~= nne
        error('real_HO_flux:dimensionMismatch', ...
            'Expected u and uPN to have %d columns (nec+1).', nne);
    end

    % M1-Fluss mit einheitlicher Wellengeschwindigkeit begrenzen.
    % PN liefert nur die HO-Momente, nicht die Charakteristikgeschwindigkeit.
    lambda = 1;

    if ilim ~= -1
        error('real_HO_flux:unsupportedCoupledMethod', ...
            ['Coupled M1PN currently supports only ilim=-1 (LLF). ' ...
            'Requested ilim=%d (%s).'], ilim, method_name(ilim));
    end

    pn_obs = pn_closure_observables(uPN(:, 1:nec), u(:, 1:nec));
    pn_m2 = reshape(pn_obs.m2, 1, []);
    pn_m2 = [pn_m2, pn_m2(1)];

    Gstar = LLF(u, dt, dx, flux, lambda, pn_m2, nec);
end

function name = method_name(ilim)
    switch ilim
        case 1
            name = 'MCL';
        case -1
            name = 'LLF';
        otherwise
            name = 'LW';
    end
end

