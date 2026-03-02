function compare_limiters(profile, keep_figures)
%COMPARE_LIMITERS Run paper-vs-MCL limiter comparison on plane-source benchmark.
%
% compare_limiters()              -> quick profile
% compare_limiters('extended')    -> broader profile
% compare_limiters('quick', true) -> keep figures open

if nargin < 1 || isempty(profile)
    profile = 'quick';
end
if nargin < 2
    keep_figures = false;
end

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'src'));

cfg = mm_default_config();
cfg.io.close_figures = ~logical(keep_figures);

switch lower(profile)
    case 'quick'
        % Keep runtime practical for direct comparison.
        cfg.models.families = {'HFMn', 'PMMn'};
        cfg.models.orders.HFMn = [4, 8];
        cfg.models.orders.PMMn = [4];

        cfg.paper2.n_cells = 40;
        cfg.paper2.tf = 0.25;

        cfg.reference.n_cells = 180;
        cfg.reference.n_mu = 64;
        cfg.reference.tf = cfg.paper2.tf;
        cfg.reference.domain = cfg.paper2.domain;

    case 'extended'
        cfg.models.families = {'PN', 'MN', 'HFPn', 'HFMn', 'PMPn', 'PMMn'};
        cfg.models.orders.PN = [3, 7];
        cfg.models.orders.MN = [3];
        cfg.models.orders.HFPn = [4, 8];
        cfg.models.orders.HFMn = [4, 8];
        cfg.models.orders.PMPn = [4, 8];
        cfg.models.orders.PMMn = [4];

        cfg.paper2.n_cells = 48;
        cfg.paper2.tf = 0.3;

        cfg.reference.n_cells = 220;
        cfg.reference.n_mu = 64;
        cfg.reference.tf = cfg.paper2.tf;
        cfg.reference.domain = cfg.paper2.domain;

        cfg.paths.results = fullfile(cfg.paths.project_root, 'results', 'compare_extended_run');

    otherwise
        error('Unknown profile ''%s''. Use ''quick'' or ''extended''.', profile);
end

cmp = compare_limiters_mcl_vs_paper(cfg);
disp(cmp.table);

end
