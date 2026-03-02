function results = runtests_mm()
%RUNTESTS_MM Run MATLAB unit tests for moment-model project.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'tests'));

results = runtests('TestMomentModels');
disp(results);

end
