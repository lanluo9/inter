% icell = 1
% icond = 3 % isi = 250 ms
% 
% data % resp to 8 ori 
% ori_pref_cond(icell) % pref ori via curve fit
% 
% 
% plot([0:22.5:180], [data, data(1)])
% hold on
% xline(ori_pref_cond(icell))

load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1323_200720_003\fit_tuning_isi3.mat')
load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1323_200720_003\fit_bootstrap_90perc.mat')

fit_param_cell = squeeze(fit_param(64, :, :));

R2_noad = fit_param(:, 7, 1); % ncell, final param, cond1 (noad)
R2_250  = fit_param(:, 7, 3); % cond3 (250 isi)

plot(R2_noad(well_fit_cell)) 
hold on
plot(R2_250(well_fit_cell))

%%

load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\LM_i1328_201119_003\fit_tuning_isi3.mat')
load('Z:\All_Staff\home\lan\Data\2P_images\mat_inter\LM_i1328_201119_003\fit_bootstrap_90perc.mat')

fit_param_cell = squeeze(fit_param(71, :, :));

R2_noad = fit_param(:, 7, 1); % ncell, final param, cond1 (noad)
R2_750  = fit_param(:, 7, 2); % cond2 (750 isi)
R2_250  = fit_param(:, 7, 3); % cond3 (250 isi)

plot(R2_noad(well_fit_cell)) 
hold on
plot(R2_750(well_fit_cell))
plot(R2_250(well_fit_cell))
legend