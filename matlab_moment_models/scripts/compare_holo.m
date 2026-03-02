function compare_holo(profile, keep_figures)
%COMPARE_HOLO Compare HOLOM1PN against baseline models on plane-source test.
%
% compare_holo()               -> quick profile
% compare_holo('extended')     -> broader profile
% compare_holo('quick', true)  -> keep figures open

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
cfg.holo.enabled = true;
cfg.holo.pn_orders = [3,7];

switch lower(profile)
    case 'quick'
        cfg.models.families = {'PN', 'MN'};
        cfg.models.orders.PN = [3];
        cfg.models.orders.MN = [1];

        cfg.paper2.n_cells = 20;
        cfg.paper2.tf = 1;

        cfg.reference.n_cells = 100;
        cfg.reference.n_mu = 40;
        cfg.reference.tf = cfg.paper2.tf;
        cfg.reference.domain = cfg.paper2.domain;

    case 'extended'
        cfg.models.families = {'PN', 'MN', 'HFPn', 'HFMn', 'PMPn', 'PMMn'};
        cfg.models.orders.PN = [3, 7];
        cfg.models.orders.MN = [1, 3];
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

        cfg.paths.results = fullfile(cfg.paths.project_root, 'results', 'compare_holo_extended_run');

    otherwise
        error('Unknown profile ''%s''. Use ''quick'' or ''extended''.', profile);
end

cmp = compare_holo_m1pn_vs_all(cfg);
disp(cmp.table);

end
