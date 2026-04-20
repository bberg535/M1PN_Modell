function ctx = mm_parallel_context(par_cfg, work_size, work_kind)
%MM_PARALLEL_CONTEXT Resolve parallel execution context and optional pool startup.
% The mode controls where parallelism is allowed:
% - single_run: inner solver loops (cells/interfaces/directions)
% - ensemble: outer independent model/case loops

if nargin < 1 || isempty(par_cfg)
    par_cfg = struct();
end
if nargin < 2 || isempty(work_size)
    work_size = 0;
end
if nargin < 3 || isempty(work_kind)
    work_kind = 'cells';
end

ctx = struct();
ctx.enabled = logical(get_field_or(par_cfg, 'enabled', false));
ctx.mode = normalize_mode(get_field_or(par_cfg, 'mode', 'single_run'));
ctx.num_workers = get_field_or(par_cfg, 'num_workers', 0);
ctx.pool_type = get_field_or(par_cfg, 'pool_type', 'local');
ctx.allow_fp_drift = logical(get_field_or(par_cfg, 'allow_fp_drift', true));
ctx.min_cells = get_field_or(par_cfg, 'min_cells', 128);
ctx.min_interfaces = get_field_or(par_cfg, 'min_interfaces', 128);
ctx.work_kind = lower(char(work_kind));
ctx.work_size = work_size;
ctx.threshold = threshold_for_kind(ctx, ctx.work_kind);
ctx.pool_started = false;
ctx.pool_available = false;
ctx.use_parallel = false;

if ~ctx.enabled
    return;
end

if ~has_parallel_support()
    return;
end

if ~ctx.allow_fp_drift && ~strcmp(ctx.work_kind, 'ensemble')
    return;
end

if ~mode_allows_kind(ctx.mode, ctx.work_kind)
    return;
end

if ctx.work_size < ctx.threshold
    return;
end

pool = gcp('nocreate');
if isempty(pool)
    pool = try_start_pool(ctx.pool_type, ctx.num_workers);
    ctx.pool_started = ~isempty(pool);
end

ctx.pool_available = ~isempty(pool);
ctx.use_parallel = ctx.pool_available;

end

function tf = mode_allows_kind(mode, work_kind)
switch mode
    case 'single_run'
        tf = strcmp(work_kind, 'cells') || strcmp(work_kind, 'interfaces') || strcmp(work_kind, 'directions');
    case 'ensemble'
        tf = strcmp(work_kind, 'ensemble');
    otherwise
        tf = false;
end
end

function thr = threshold_for_kind(ctx, work_kind)
switch work_kind
    case 'cells'
        thr = ctx.min_cells;
    case 'interfaces'
        thr = ctx.min_interfaces;
    case 'directions'
        thr = ctx.min_cells;
    case 'ensemble'
        thr = 2;
    otherwise
        thr = inf;
end
end

function pool = try_start_pool(pool_type, num_workers)
persistent start_failed;

if isempty(start_failed)
    start_failed = false;
end

pool = gcp('nocreate');
if ~isempty(pool) || start_failed
    return;
end

try
    if num_workers > 0
        pool = parpool(pool_type, num_workers);
    else
        pool = parpool(pool_type);
    end
catch
    start_failed = true;
    pool = [];
end
end

function tf = has_parallel_support()
tf = (exist('gcp', 'file') == 2 || exist('gcp', 'builtin') == 5) && ...
     (exist('parpool', 'file') == 2 || exist('parpool', 'builtin') == 5);
end

function mode = normalize_mode(in)
if isstring(in)
    mode = lower(char(in));
elseif ischar(in)
    mode = lower(strtrim(in));
else
    mode = 'single_run';
end

if ~strcmp(mode, 'single_run') && ~strcmp(mode, 'ensemble')
    mode = 'single_run';
end
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
