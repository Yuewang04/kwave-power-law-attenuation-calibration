function run_all_benchmarks(results_dir)
%RUN_ALL_BENCHMARKS Single-session master runner for all 1-D k-Wave benchmarks.
%   Executes every 1-D experiment in one MATLAB process to avoid repeated JVM
%   startup overhead and regenerates every 1-D JSON summary in results/.

if nargin < 1 || isempty(results_dir)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    results_dir = fullfile(root, 'results');
end
if ~isfolder(results_dir), mkdir(results_dir); end

t_total = tic;
fprintf('============================================================\n');
fprintf('Starting Complete 1-D k-Wave Benchmark Suite\n');
fprintf('Output directory: %s\n', results_dir);
fprintf('============================================================\n\n');

% 1. Primary y=2 experiment -> matlab_results.json and results.json
fprintf('=== [A] Primary y=2 Calibration & Frequency Sweep ===\n');
raw_path = fullfile(results_dir, 'matlab_results.json');
res_y2 = run_kwave_experiments(raw_path, 2);
copyfile(raw_path, fullfile(results_dir, 'results.json'));

% 2. Exponent sweep -> results_n_sweep.json
fprintf('\n=== [B] Exponent Sweep (y = 0.5 ... 3.0) ===\n');
exponents = [0.5, 0.9, 1.0, 1.05, 1.1, 1.5, 2.0, 2.1, 2.5, 3.0];
sweep_cells = cell(1, numel(exponents));
for idx_e = 1:numel(exponents)
    ev = exponents(idx_e);
    fprintf('  -> Exponent y = %.2f\n', ev);
    tmp_path = fullfile(results_dir, sprintf('_tmp_n_%s.json', strrep(num2str(ev), '.', 'p')));
    try
        if ev == 2.0
            one_res = res_y2;
        else
            one_res = run_kwave_experiments(tmp_path, ev);
        end
        if isfile(tmp_path), delete(tmp_path); end
        one_res.n = ev;
        sweep_cells{idx_e} = one_res;
    catch ME
        if isfile(tmp_path), delete(tmp_path); end
        sweep_cells{idx_e} = struct( ...
            'n', ev, ...
            'status', 'unsupported_by_kwave', ...
            'error_identifier', ME.identifier, ...
            'error_message', ME.message, ...
            'reason', sprintf('k-Wave input validation rejected alpha_power=%.4g (%s). Note tan(pi*y/2) also has poles at odd integers y=1, 3.', ev, ME.message));
    end
end
% Encode cell array directly so heterogeneous structs (valid runs + y=3 rejection) serialize as JSON array
fid = fopen(fullfile(results_dir, 'results_n_sweep.json'), 'w');
assert(fid >= 0, 'Cannot open results_n_sweep.json');
fwrite(fid, jsonencode(sweep_cells, 'PrettyPrint', true), 'char');
fclose(fid);

% 3. Near-y=1 multi-frequency scans (y = 0.98, 0.99, 1.01, 1.02)
fprintf('\n=== [C] Near-y=1 Multi-Frequency Scans ===\n');
run_near_y_frequency_scans(results_dir);

% 4. Timestep & estimator controls (y = 0.98 and y = 0.99)
fprintf('\n=== [D] Timestep & Estimator Controls ===\n');
run_controls(fullfile(results_dir, 'controls_y098.json'));

% 5. Fine singularity scan, alpha_coeff scaling, grid convergence, Part 3, fitPowerLawParams audit
fprintf('\n=== [E] Singularity Scan, Scaling, Grid, Worked Example & fitPowerLawParams Audit ===\n');
run_singularity_and_diagnostics(results_dir);

% 6. Also refresh interior_source_validation.json for backward compatibility
validate_interior_source(fullfile(results_dir, 'interior_source_validation.json'));

fprintf('\n============================================================\n');
fprintf('ALL 1-D BENCHMARKS COMPLETED IN %.1f s!\n', toc(t_total));
fprintf('============================================================\n');
end
