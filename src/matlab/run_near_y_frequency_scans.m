function run_near_y_frequency_scans(output_dir)
% Extend the original protocol to near-y=1 exponents across center frequencies.
% Coefficients remain matched at 0.5 MHz. This is an exploratory baseline.
if nargin < 1
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
