function run_near_y_frequency_scans(output_dir)
%RUN_NEAR_Y_FREQUENCY_SCANS Near-y=1 frequency sweeps with standard interior source.
%   Generates results/near_y_frequency_y{0p98,0p99,1p01,1p02}.json.

if nargin < 1 || isempty(output_dir)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    output_dir = fullfile(root, 'results');
end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

for y = [0.98 0.99 1.01 1.02]
    tag = strrep(sprintf('%.2f', y), '.', 'p');
    output = fullfile(output_dir, ['near_y_frequency_y' tag '.json']);
    run_kwave_experiments(output, y);
end
end
