function [tuning_curve_cell_isi] = tuning_curve_nofit(dfof_tg, save_flag)

% save tuning curve data for cells, using resp to target
% input: dfof_tg = ncell x nori x nisi [noad ad750 ad250]. save_flag to toggle save .mat
% output: 
% tuning_curve_cell_isi = ncell x nstim x nisi [noad vs ad750 vs ad250]

global ncell ori_list nori

ncond = size(dfof_tg, 3);
tuning_curve_cell_isi = pi * ones(ncell, ncond, nori);
    
for icond = 1:ncond
    dfof_cond = dfof_tg(:,:,icond); 
    for icell = 1 : ncell
        ori_rad = deg2rad(ori_list);
        data = dfof_cond(icell,:);
        tuning_curve_cell_isi(icell, icond, :) = data; % tuning curve across 8 ori
    end
end

if save_flag; save tuning_curve.mat tuning_curve_cell_isi; end

