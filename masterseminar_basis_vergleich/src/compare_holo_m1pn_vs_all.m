function cmp = compare_holo_m1pn_vs_all(cfg)
%COMPARE_HOLO_M1PN_VS_ALL Compare coupled HOLOM1PN against baseline models.

if nargin < 1 || isempty(cfg)
    cfg = mm_default_config();
end

if ~exist(cfg.paths.results, 'dir')
    mkdir(cfg.paths.results);
end

closeFigures = get_close_figures(cfg);
includeBaselines = logical(get_field_or(cfg.holo, 'compare_include_baselines', true));

resBase = struct('table', table());
if includeBaselines
    cfgBase = cfg;
    cfgBase.limiter.type = 'mcl';
    cfgBase.models = ensure_m1_baseline(cfgBase.models);
    resBase = run_paper2_plane_source(cfgBase);
end

resHolo = run_holo_m1pn_plane_source(cfg);

baseVars = {'model', 'order', 'nMom', 'L1', 'Linf', 'symmetry_L1_rel', 'min_rho', 'runtime_s', 'n_steps'};
Tb = ensure_columns(resBase.table, baseVars);
Th = ensure_columns(resHolo.table, baseVars);

if isempty(Tb)
    Tall = Th;
else
    Tall = [Tb; Th];
end

outCsv = fullfile(cfg.paths.results, 'holo_m1pn_vs_all.csv');
writetable(Tall, outCsv);

plot_comparison(Tall, cfg.paths.results, closeFigures);

cmp = struct();
cmp.table = Tall;
cmp.csv = outCsv;
cmp.baselines = resBase;
cmp.holo = resHolo;

end

function models = ensure_m1_baseline(models)
if ~isfield(models.orders, 'MN') || isempty(models.orders.MN)
    models.orders.MN = 1;
end

mnOrders = models.orders.MN(:).';
if ~ismember(1, mnOrders)
    models.orders.MN = [1, mnOrders];
end

if ~ismember('MN', models.families)
    models.families = [{'MN'}, models.families];
end
end

function T = ensure_columns(Tin, vars)
if isempty(Tin)
    T = table();
    for i = 1:numel(vars)
        name = vars{i};
        if strcmp(name, 'model')
            T.(name) = strings(0, 1);
        else
            T.(name) = zeros(0, 1);
        end
    end
    return;
end

T = Tin;
for i = 1:numel(vars)
    name = vars{i};
    if ~ismember(name, T.Properties.VariableNames)
        if strcmp(name, 'model')
            T.(name) = repmat("", height(T), 1);
        else
            T.(name) = nan(height(T), 1);
        end
    end
end
T = T(:, vars);
end

function plot_comparison(T, outDir, closeFigures)
labels = strcat(string(T.model), "-", string(T.order));

fig = figure('Color', 'w');
tiledlayout(1, 3);

nexttile;
bar(categorical(labels), T.L1);
ylabel('L1 error');
title('L1 comparison');
grid(gca, 'on');

nexttile;
bar(categorical(labels), T.Linf);
ylabel('Linf error');
title('Linf comparison');
grid(gca, 'on');

nexttile;
bar(categorical(labels), T.runtime_s);
ylabel('runtime [s]');
title('Runtime comparison');
grid(gca, 'on');

outBase = fullfile(outDir, 'holo_m1pn_vs_all');
saveas(fig, [outBase '.png']);
try
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
catch
end
if closeFigures
    close(fig);
end
end

function closeFigures = get_close_figures(cfg)
closeFigures = true;
if isfield(cfg, 'io') && isfield(cfg.io, 'close_figures')
    closeFigures = logical(cfg.io.close_figures);
end
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
