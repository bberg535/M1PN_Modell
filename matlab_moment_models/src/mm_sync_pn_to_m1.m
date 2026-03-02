function [u_pn_sync, sync_state] = mm_sync_pn_to_m1(u_pn_ho, u_m1_target, sync_cfg)
%MM_SYNC_PN_TO_M1 Project PN state onto constraint u(1:2,:)=u_m1_target.

if nargin < 3 || isempty(sync_cfg)
    sync_cfg = struct();
end
tSync = tic;

mode = lower(get_field_or(sync_cfg, 'weight', 'identity'));

[nMom, nCells] = size(u_pn_ho);
if nMom < 2
    error('PN state must contain at least two moments.');
end
if size(u_m1_target, 1) ~= 2 || size(u_m1_target, 2) ~= nCells
    error('u_m1_target must be of size [2, nCells].');
end

switch mode
    case 'identity'
        % Closed-form constrained least-squares projection in Euclidean norm.
        u_pn_sync = u_pn_ho;
        u_pn_sync(1:2, :) = u_m1_target;
    otherwise
        error('Unsupported sync weight: %s', mode);
end

delta = u_pn_sync - u_pn_ho;
constraint_violation = max(abs(u_pn_sync(1:2, :) - u_m1_target), [], 'all');

sync_state = struct();
sync_state.mode = mode;
sync_state.correction_norm_fro = norm(delta, 'fro');
sync_state.correction_norm_max = max(abs(delta), [], 'all');
sync_state.max_constraint_violation = constraint_violation;
sync_state.elapsed_s = toc(tSync);

end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
