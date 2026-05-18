function m2 = m1pn_second_moment(uPN)
%M1PN_SECOND_MOMENT Return the PN-based second moment used by coupled M1PN.

    obs = pn_observables(uPN);
    m2 = reshape(obs.m2, 1, []);
end
