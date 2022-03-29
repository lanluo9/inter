function well_fit_cell = well_fit_cell_criteria(dfof_align_tg, nrun, save_flag)

% define well-fit: noad-tg, 90% bootstraps within 22.5 deg of all-trials-included fit
% default bootstrap nrun = 1000
% input: dfof_align_tg. output: well_fit_cell as logical

global ncell nori id_ori ori_list id_noad range_base range_resp
dfof_avg_runs = pi * ones(ncell, nori, nrun); % dfof_ste_runs = pi * ones(ncell, nori, nrun);
fit_param_runs = pi * ones(ncell, 7, nrun); ori_pref_runs = pi * ones(ncell, nrun);

if nargin == 1; nrun = 1000; end
disp('start bootstrap runs')
for irun = 1 : nrun
    if ~mod(irun, 100); disp(num2str(irun)); end

for icell = 1 : ncell        
    for iori = 1 : nori
        idx = intersect(id_ori{iori}, id_noad);
        ntrials_ori_noad = length(idx);
        bootstrap_draw = round(ntrials_ori_noad * 0.7);
        idx_run = randsample(idx, bootstrap_draw, 1); % w replacement

        base_win = squeeze(dfof_align_tg(icell, idx_run, range_base)); base_win = mean(base_win, 2);
        resp_win = squeeze(dfof_align_tg(icell, idx_run, range_resp)); resp_win = mean(resp_win, 2);
        dfof_avg_runs(icell, iori, irun) = mean( resp_win - base_win );
%       dfof_ste_runs(icell, iori, irun) = std( resp_win - base_win ) ./ sqrt(ntrials_ori_noad);
    end

    data = dfof_avg_runs(icell, :, irun); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(deg2rad(ori_list), data);
    fit_param_runs(icell, :, irun) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
%   icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2

    ori_pref = rad2deg(u1_hat);
    ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180; ori_pref(ori_pref >= 180) = ori_pref(ori_pref >= 180) - 180;
    ori_pref_runs(icell, irun) = ori_pref;
end
end

ori_distance = abs(ori_pref_runs - mean(ori_pref_runs,2));
ori_distance(ori_distance > 90) = ori_distance(ori_distance > 90) - 90; % distance btw oris <= 90. flip those too large
ori_distance = sort(ori_distance, 2); % sort each row as cell

percentile_threshold = 0.90; percentile_idx = percentile_threshold * nrun;
ori_perc = ori_distance(:, percentile_idx);
well_fit_cell = ori_perc < (180 / nori);

if save_flag
    save fit_bootstrap_90perc.mat fit_param_runs ori_pref_runs well_fit_cell ori_perc; 
end


