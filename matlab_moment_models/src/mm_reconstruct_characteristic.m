function [uL, uR, rec_state] = mm_reconstruct_characteristic(u_cell, jac_state, grid, rec_cfg)
%MM_RECONSTRUCT_CHARACTERISTIC Piecewise-linear reconstruction in (optional) characteristic vars.
% Reference: Seminarquelle 2 Eq. (4.16) (minmod reconstruction, preferably in characteristic variables).

u = u_cell;
[nMom, nCells] = size(u);
dz = grid.dz;

if isfield(rec_cfg, 'use_characteristic')
    useChar = logical(rec_cfg.use_characteristic);
else
    useChar = true;
end

slopes = zeros(nMom, nCells);

for i = 2:(nCells - 1)
    duF = u(:, i + 1) - u(:, i);
    duB = u(:, i) - u(:, i - 1);
    duC = 0.5 * (u(:, i + 1) - u(:, i - 1));

    if useChar && ~isempty(jac_state) && isfield(jac_state, 'V') && numel(jac_state.V) >= i
        V = jac_state.V{i};
        Vinv = jac_state.Vinv{i};
    else
        V = eye(nMom);
        Vinv = eye(nMom);
    end

    cF = Vinv * duF;
    cB = Vinv * duB;
    cC = Vinv * duC;

    cSlope = mm_minmod(cF, cB, cC) / dz;
    slopes(:, i) = V * cSlope;
end

uL = u - 0.5 * dz * slopes;
uR = u + 0.5 * dz * slopes;

% First-order at boundaries.
uL(:, 1) = u(:, 1);
uR(:, 1) = u(:, 1);
uL(:, nCells) = u(:, nCells);
uR(:, nCells) = u(:, nCells);

rec_state = struct('slopes', slopes, 'used_characteristic', useChar);

end
