function [R_sq, ori_fit] = fit_tuning_jeff(dfof_tg, save_flag)

% write data in correct format for jeff population vector decoder - jin2019
% input: 
    % dfof_tg = ncell x nori. save_flag to toggle save .mat
% output: 
    % R_sq = ncell x 1. von mises goodness of fit
    % ori_fit = 181 deg x ncell. von mises fit tuning curve values across deg

global ncell ori_list

fit_param = pi * ones(ncell, 7);
ori_fit = pi * ones(181, ncell);
dfof_cond = dfof_tg(:,:,1); % jeff code only fits unadapted resp to von mises

theta_deg_arr = (0:1:180)';
theta_rad_arr = deg2rad(theta_deg_arr);

for icell = 1 : ncell
    ori_rad = deg2rad(ori_list);
    data = dfof_cond(icell,:);

    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_squared] = miaovonmisesfit_ori(ori_rad, data);
    fit_param(icell, :) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_squared];
    % icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, 
    % sse sum of squared error, R2

    y_pred = b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(theta_rad_arr - u1_hat)) - 1)); 
    % copied from miaovonmisesfit_ori.m
    ori_fit(:, icell) = y_pred;
end
R_sq = fit_param(:, end); % aka R_squared

if save_flag; save fit_tuning_jeff.mat R_sq ori_fit; end

