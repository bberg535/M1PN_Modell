function method = select_flux_method(ilim)
%SELECT_FLUX_METHOD Map ilim to the requested transport flux.

    switch ilim
        case 1
            method = @MCL;
        case -1
            method = @LLF;
        otherwise
            method = @LW;
    end
end
