function out = mm_plot_seminarquelle2_profiles(results_dir, panels, out_base, io_cfg)
%MM_PLOT_SEMINARQUELLE2_PROFILES Create paper-like profile figures for the plane-source benchmark.

if nargin < 4 || isempty(io_cfg)
    io_cfg = struct();
end

fig = create_figure(io_cfg);
if isempty(fig)
    error('Could not create figure for Seminarquelle 2 profile plot.');
end
set_figure_size(fig, numel(panels));

tl = tiledlayout(fig, 1, numel(panels), 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(panels)
    ax = nexttile(tl);
    plot_single_panel(ax, results_dir, panels(i));
end

write_figure_files(fig, out_base, io_cfg);
if get_field_or(io_cfg, 'close_figures', true)
    close(fig);
end

out = struct();
out.base = out_base;
out.png = [out_base '.png'];
if get_field_or(io_cfg, 'write_pdf', true)
    out.pdf = [out_base '.pdf'];
else
    out.pdf = '';
end

end

function plot_single_panel(ax, results_dir, panel)
hold(ax, 'on');

orders = panel.orders;
curve_handles = gobjects(1, numel(orders));
legend_labels = cell(1, numel(orders) + 1);
ref_handle = gobjects(1, 1);

for k = 1:numel(orders)
    prof = read_profile(results_dir, panel.family, orders(k));
    mask = prof.z >= -1.0e-12;
    z = prof.z(mask);
    rho = prof.rho_model(mask);
    rho_ref = prof.rho_ref(mask);

    if k == 1
        ref_handle = plot(ax, z, rho_ref, 'k--', 'LineWidth', 1.1, 'DisplayName', 'Ref');
    end

    sty = order_style(k);
    marker_idx = unique(max(1, round(linspace(1, numel(z), min(6, numel(z))))));
    curve_handles(k) = plot(ax, z, rho, ...
        'Color', sty.color, ...
        'LineStyle', '-', ...
        'LineWidth', 1.1, ...
        'Marker', sty.marker, ...
        'MarkerSize', 4.8, ...
        'MarkerIndices', marker_idx, ...
        'DisplayName', model_curve_label(panel.family, orders(k)));
    legend_labels{k} = model_curve_label(panel.family, orders(k));
end

legend_labels{end} = 'Ref';
legend(ax, [curve_handles, ref_handle], legend_labels, ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal', ...
    'Interpreter', 'tex');

grid(ax, 'on');
box(ax, 'on');
ax.LineWidth = 0.6;
ax.TickDir = 'in';
ax.FontSize = 11;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

xlim(ax, [0, 1]);
xticks(ax, 0:0.2:1);
ylim(ax, [-1.2, 12]);
yticks(ax, 0:2:12);

xlabel(ax, panel.x_label, 'Interpreter', 'tex');
ylabel(ax, panel.y_label, 'Interpreter', 'tex');
add_panel_title(ax, panel.panel_title);

end

function T = read_profile(results_dir, model_name, order)
fn = fullfile(results_dir, sprintf('paper2_plane_source_%s_%d_profile.csv', model_name, order));
if ~exist(fn, 'file')
    error('Missing profile file: %s', fn);
end
T = readtable(fn);
end

function sty = order_style(idx)
styles = { ...
    struct('color', [0.2, 0.45, 0.9], 'marker', 'o'), ...
    struct('color', [0.95, 0.65, 0.05], 'marker', 'o'), ...
    struct('color', [0.0, 0.65, 0.2], 'marker', 'd')};
idx = max(1, min(idx, numel(styles)));
sty = styles{idx};
end

function s = model_curve_label(model_name, order)
switch char(model_name)
    case 'PN'
        prefix = 'P';
    case 'HFMn'
        prefix = 'HFM';
    case 'HFPn'
        prefix = 'HFP';
    case 'PMPn'
        prefix = 'PMP';
    otherwise
        prefix = char(model_name);
end
s = sprintf('%s%d', prefix, order);
end

function fig = create_figure(io_cfg)
visible_state = ternary(get_field_or(io_cfg, 'close_figures', true), 'off', 'on');
fig = [];
try
    fig = figure('Color', 'w', 'Visible', visible_state);
catch
    try
        set(groot, 'defaultFigureVisible', visible_state);
        fig = figure('Color', 'w');
        set(fig, 'Visible', visible_state);
    catch
        fig = [];
    end
end
end

function set_figure_size(fig, n_panels)
if n_panels <= 1
    fig.Position(3:4) = [760, 520];
else
    fig.Position(3:4) = [1180, 520];
end
end

function add_panel_title(ax, panel_title)
text(ax, 0.5, 0.99, panel_title, ...
    'Units', 'normalized', ...
    'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12, ...
    'Clipping', 'off');
end

function write_figure_files(fig, out_base, io_cfg)
try
    exportgraphics(fig, [out_base '.png'], 'Resolution', 200);
catch
    try
        print(fig, [out_base '.png'], '-dpng', '-r200');
    catch
    end
end

if get_field_or(io_cfg, 'write_pdf', true)
    try
        exportgraphics(fig, [out_base '.pdf'], 'ContentType', 'vector');
    catch
        try
            print(fig, [out_base '.pdf'], '-dpdf', '-painters');
        catch
        end
    end
end
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function out = ternary(cond, true_value, false_value)
if cond
    out = true_value;
else
    out = false_value;
end
end
