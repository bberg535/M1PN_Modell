function run_model()
%RUN_ALL Execute both benchmark scripts with default configuration.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'src'));

% Ensure workers do not keep stale function versions across reruns.
if exist('gcp', 'file') == 2 || exist('gcp', 'builtin') == 5
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
end

cfg = mm_default_config();
cfg.models.families = {'PMMn'}
cfg.models.orders.PMMn = [4, 8];
cfg.parallel.enabled = true;
cfg.parallel.mode = 'single-run';
cfg.parallel.pool_type = 'local';
cfg.parallel.num_workers = 8;
cfg.parallel.force_serial_partial_entropy = true;

% Make PMMn behavior explicit: component-wise reconstruction is currently
% more robust on the plane-source benchmark than characteristic PMMn.
cfg.reconstruction.use_characteristic_partial_entropy = false;
cfg.optimizer.use_change_of_basis_partial_entropy = false;

fprintf('[run_model] use_characteristic_partial_entropy = %d\n', cfg.reconstruction.use_characteristic_partial_entropy);
fprintf('[run_model] use_change_of_basis_partial_entropy = %d\n', cfg.optimizer.use_change_of_basis_partial_entropy);
fprintf('[run_model] force_serial_partial_entropy = %d\n', cfg.parallel.force_serial_partial_entropy);
fprintf('[run_model] cfl_safety = %.3f\n', cfg.solver.cfl_safety);
if isfield(cfg.solver, 'cfl_safety_partial_entropy')
    fprintf('[run_model] cfl_safety_partial_entropy = %.3f\n', cfg.solver.cfl_safety_partial_entropy);
end
fprintf('[run_model] mm_step_flux_rk2 from: %s\n', which('mm_step_flux_rk2'));
fprintf('[run_model] mm_default_config from: %s\n', which('mm_default_config'));


%disp('Running paper 1 (Gauss/Heaviside)...');
%r1 = run_paper1_gauss_heaviside(cfg);
%disp(r1.table);

disp('Running paper 2 (Plane source)...');
r2 = run_paper2_plane_source(cfg);
disp(r2.table);

end
