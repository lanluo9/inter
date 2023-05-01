function [fit_param, ori_pref, tuning_curve_cell_isi] = fit_tuning(dfof_tg, save_flag)

% find tuning curve fit & preferred orientation for cells, using resp to target
% input: dfof_tg = ncell x nori x nisi [noad ad750 ad250]. save_flag to toggle save .mat
% output: 
% fit_param = ncell x nparam x nisi [noad vs ad750 vs ad250]
% ori_pref = ncell x nisi [noad vs ad750 vs ad250]

global ncell ori_list nori

ncond = size(dfof_tg, 3);
fit_param = pi * ones(ncell, 7, ncond);
ori_pref = pi * ones(ncell, ncond);
tuning_curve_cell_isi = pi * ones(ncell, ncond, nori);
    
for icond = 1:ncond
    dfof_cond = dfof_tg(:,:,icond); 
    for icell = 1 : ncell
        ori_rad = deg2rad(ori_list);
        data = dfof_cond(icell,:);
        tuning_curve_cell_isi(icell, icond, :) = data; % tuning curve across 8 ori

        [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(ori_rad, data);
        fit_param(icell, :, icond) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
        % icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, 
        % sse sum of squared error, R2
    end
    ori_pref_cond = rad2deg(fit_param(:, 5, icond)); % = u1_hat
    ori_pref_cond(ori_pref_cond < 0) = ori_pref_cond(ori_pref_cond < 0) + 180; 
    ori_pref_cond(ori_pref_cond > 180) = ori_pref_cond(ori_pref_cond > 180) - 180;
    ori_pref(:, icond) = ori_pref_cond;
end

if save_flag; save fit_tuning_isi3.mat fit_param ori_pref tuning_curve_cell_isi; end

