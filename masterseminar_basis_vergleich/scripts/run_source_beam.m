function res = run_source_beam(mode)
%RUN_SOURCE_BEAM Execute source-beam benchmark (Seminarquelle 2, Section 6.1.2).
% mode:
%   'quick'   - default model orders from mm_default_config.
%   'figure5' - convergence sweeps comparable to Figure 5.

if nargin < 1 || isempty(mode)
    mode = 'quick';
end

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
cfg.models.families = {'HFMn', 'HFPn', 'PMMn', 'PMPn', 'MN', 'PN'};

% Preserve previous script behavior: run with local process pool by default.
cfg.parallel.enabled = true;
cfg.parallel.mode = 'single_run';
cfg.parallel.pool_type = 'local';
cfg.parallel.num_workers = 8;

switch lower(char(mode))
    case 'figure5'
        cfg.models.orders.MN = 3;
        cfg.models.orders.PN = [3,7];
        cfg.models.orders.HFMn = 4;
        cfg.models.orders.HFPn = [4,8];
        cfg.models.orders.PMMn = 4;
        cfg.models.orders.PMPn = [4,8];
    case 'quick'
        % Keep default orders.
    otherwise
        error('Unknown mode "%s". Use "quick" or "figure5".', mode);
end

disp(['Running paper 2 (Source beam), mode=' char(mode) ' ...']);
res = run_paper2_source_beam(cfg);

if istable(res.table)
    keepCols = {'model', 'order', 'nMom', 'L1', 'Linf'};
    keepCols = keepCols(ismember(keepCols, res.table.Properties.VariableNames));
    disp(res.table(:, keepCols));
else
    disp(res.table);
end

end
