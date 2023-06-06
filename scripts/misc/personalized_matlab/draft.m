
tmp = squeeze(nanmean(nanmean(dfof_align_ad, 1), 2));
% tmp = squeeze(dfof_align_ad(25, 16, :));
figure
plot(tmp)
hold on
yline(0)
% xlim([0, 60]);

% 10(11), 20(19), 34
% 8 or 10 frame = 250 ms
% 24 frame = 750 ms

%%

figure
imagesc(tc_trial_base_avg)

%%


% % % icell = 1
% % % icond = 3 % isi = 250 ms
% % % 
% % % data % resp to 8 ori 
% % % ori_pref_cond(icell) % pref ori via curve fit
% % % 
% % % 
% % % plot([0:22.5:180], [data, data(1)])
% % % hold on
% % % xline(ori_pref_cond(icell))
% % 
% % % load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1323_200720_003\fit_tuning_isi3.mat')
% % % load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1323_200720_003\fit_bootstrap_90perc.mat')
% % % load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1323_200720_003\dfof_trial.mat')
% % 
% % cell_id = 64 % 1-based indexing
% % 
% % fit_param_cell = squeeze(fit_param(cell_id, :, :));
% % 
% % R2_noad = fit_param(:, 7, 1); % ncell, final param, cond1 (noad)
% % R2_250  = fit_param(:, 7, 3); % cond3 (250 isi)
% % 
% % figure;
% % plot(R2_noad(well_fit_cell)) 
% % hold on
% % plot(R2_250(well_fit_cell))
% % 
% % %%
% % % resp_cell = dfof_ad_trial(cell_id, :, :) - dfof_base_trial(cell_id, :, :);
% % resp_cell_noad = resp_cell(:, 1);
% % resp_cell_250 = resp_cell(:, 3);
% % 
% % figure;
% % plot(resp_cell_noad) 
% % hold on
% % plot(resp_cell_250)
% % 
% % %%
% % 
% % % load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\LM_i1328_201119_003\fit_tuning_isi3.mat')
% % % load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\LM_i1328_201119_003\fit_bootstrap_90perc.mat')
% % % 
% % % fit_param_cell = squeeze(fit_param(71, :, :));
% % % 
% % % R2_noad = fit_param(:, 7, 1); % ncell, final param, cond1 (noad)
% % % R2_750  = fit_param(:, 7, 2); % cond2 (750 isi)
% % % R2_250  = fit_param(:, 7, 3); % cond3 (250 isi)
% % % 
% % % plot(R2_noad(well_fit_cell)) 
% % % hold on
% % % plot(R2_750(well_fit_cell))
% % % plot(R2_250(well_fit_cell))
% % % legend
% 
% %%
% icell = 64;
% % icell = 71;
% 
% hist(fit_param(:, end, end))
% 
% figure;
% subplot(1,2,1)
% 
% icond = 1;
% dfof_cond = dfof_tg(:,:,icond); 
% ori_rad = deg2rad(ori_list);
% data = dfof_cond(icell,:); 
% [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(ori_rad, data);
% 
% theta_finer = deg2rad(0:1:179);
% y_fit = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(theta_finer - u1_hat))-1));
% 
% plot(ori_list, data)
% hold on
% plot(rad2deg(theta_finer), y_fit, '--')
% 
% subplot(1,2,2)
% 
% icond = 3;
% dfof_cond = dfof_tg(:,:,icond); 
% ori_rad = deg2rad(ori_list);
% data = dfof_cond(icell,:); 
% [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(ori_rad, data);
% 
% theta_finer = deg2rad(0:1:179);
% y_fit = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(theta_finer - u1_hat))-1));
% 
% plot(ori_list, data)
% hold on
% plot(rad2deg(theta_finer), y_fit, '--')