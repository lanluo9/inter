function well_fit_cell = well_fit_cell_criteria_relax(percentile_threshold, nrun, save_flag)

% reload ori_pref_runs from fit_bootstrap_90perc.mat 
% redefine well-fit cells as less percentile of fit pref ori falling within 22.5 deg

global nori
tmp = load("fit_bootstrap_90perc.mat", "ori_pref_runs");
ori_pref_runs = tmp.ori_pref_runs;

ori_distance = abs(ori_pref_runs - mean(ori_pref_runs,2));
ori_distance(ori_distance > 90) = ori_distance(ori_distance > 90) - 90; % distance btw oris <= 90. flip those too large
ori_distance = sort(ori_distance, 2); % sort each row as cell

percentile_idx = percentile_threshold * nrun;
ori_perc = ori_distance(:, percentile_idx);
well_fit_cell = ori_perc < (180 / nori);

if save_flag
    save fit_bootstrap_relax.mat well_fit_cell; 
end


