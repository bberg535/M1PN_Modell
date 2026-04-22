function out = mm_plot_seminarquelle2_convergence(T, families, out_base, io_cfg)
%MM_PLOT_SEMINARQUELLE2_CONVERGENCE Create paper-like convergence plot for the plane-source benchmark.

if nargin < 2 || isempty(families)
    families = {'HFMn', 'HFPn', 'PMPn', 'PN'};
end
if nargin < 4 || isempty(io_cfg)
    io_cfg = struct();
end

fig = create_figure(io_cfg);
if isempty(fig)
    error('Could not create figure for Seminarquelle 2 convergence plot.');
end
fig.Position(3:4) = [1180, 560];

tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl);
[handles, labels] = plot_metric(ax1, T, families, 'L1', 'E_h^1', '(a) L^1 convergence');
add_slope_guide(ax1, 0.2, 0.52, 0.70, '0.2th order', -22);

ax2 = nexttile(tl);
plot_metric(ax2, T, families, 'Linf', 'E_h^\infty', '(b) L^\infty convergence');
add_slope_guide(ax2, 0.8, 0.70, 0.63, '0.8th order', -28);

if ~isempty(handles)
    lgd = legend(ax1, handles, labels, ...
        'Location', 'northoutside', ...
        'Orientation', 'horizontal', ...
        'Interpreter', 'tex');
    try
        lgd.Layout.Tile = 'north';
    catch
    end
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

function [handles, labels] = plot_metric(ax, T, families, metric_name, ylabel_txt, title_txt)
hold(ax, 'on');

handles = gobjects(0);
labels = {};

for i = 1:numel(families)
    family = families{i};
    rows = T(strcmp(string(T.model), string(family)), :);
    if isempty(rows)
        continue;
    end

    [x, idx] = sort(rows.nMom);
    y = rows.(metric_name)(idx);
    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;
    x = x(valid);
    y = y(valid);
    if isempty(x)
        continue;
    end

    sty = family_style(family);
    h = loglog(ax, x, y, ...
        'Color', sty.color, ...
        'LineStyle', sty.line, ...
        'Marker', sty.marker, ...
        'LineWidth', 1.1, ...
        'MarkerSize', 5.2, ...
        'DisplayName', sty.label);

    handles(end + 1) = h; %#ok<AGROW>
    labels{end + 1} = sty.label; %#ok<AGROW>
end

grid(ax, 'on');
box(ax, 'on');
ax.LineWidth = 0.6;
ax.TickDir = 'in';
ax.FontSize = 11;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

xlabel(ax, 'n', 'Interpreter', 'tex');
ylabel(ax, ylabel_txt, 'Interpreter', 'tex');
add_panel_title(ax, title_txt);

end

function sty = family_style(family)
switch char(family)
    case 'HFMn'
        sty = struct('color', [0.2, 0.45, 0.9], 'line', '-', 'marker', 'o', 'label', 'HFM_n');
    case 'HFPn'
        sty = struct('color', [0.2, 0.45, 0.9], 'line', '--', 'marker', 'o', 'label', 'HFP_n');
    case 'PMPn'
        sty = struct('color', [0.95, 0.65, 0.05], 'line', '--', 'marker', '^', 'label', 'PMP_n');
    case 'PN'
        sty = struct('color', [0.0, 0.65, 0.2], 'line', '--', 'marker', 'p', 'label', 'P_N');
    otherwise
        sty = struct('color', [0.1, 0.1, 0.1], 'line', '-', 'marker', 'o', 'label', char(family));
end
end

function add_slope_guide(ax, slope, x_anchor_frac, y_anchor_frac, label_txt, rotation_deg)
xl = xlim(ax);
yl = ylim(ax);
if any(~isfinite([xl, yl])) || xl(1) <= 0 || yl(1) <= 0
    return;
end

logx = log(xl);
logy = log(yl);
anchor_x = exp(logx(1) + x_anchor_frac * (logx(2) - logx(1)));
anchor_y = exp(logy(1) + y_anchor_frac * (logy(2) - logy(1)));

x1 = exp(logx(1) + 0.18 * (logx(2) - logx(1)));
x2 = exp(logx(1) + 0.98 * (logx(2) - logx(1)));
x_ref = logspace(log10(x1), log10(x2), 64);
y_ref = anchor_y * (x_ref / anchor_x) .^ (-slope);

loglog(ax, x_ref, y_ref, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

label_idx = max(2, round(0.60 * numel(x_ref)));
text(ax, x_ref(label_idx), y_ref(label_idx) * 1.05, label_txt, ...
    'Interpreter', 'tex', ...
    'Rotation', rotation_deg, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);
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

function add_panel_title(ax, title_txt)
text(ax, 0.5, 1.04, title_txt, ...
    'Units', 'normalized', ...
    'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
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
