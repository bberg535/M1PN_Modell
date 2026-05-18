clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(script_dir));
setup_project_paths();

pn_debug_report('orders', [3, 5, 7], 'plane_source_method', -2, 'verbose', true);
