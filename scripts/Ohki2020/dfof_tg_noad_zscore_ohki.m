function [dfof_tg_noad, z_score] = dfof_tg_noad_zscore_ohki(dfof_align, save_flag)

% output: dfof_tg_noad & z_score = ncell x ntrial

global ncell id_isi3 range_base range_resp
base_cond = cell(ncell,1); resp_cond = cell(ncell,1);
dfof_tg_noad = cell(ncell,1); z_score_ = cell(ncell,1);

idx = id_isi3{1}; ntrial_noad = length(idx);
z_score = zeros(ncell, ntrial_noad);

for icell = 1 : ncell
    base_cond{icell,1} = mean(squeeze(dfof_align(icell, idx, range_base)),2); 
    resp_cond{icell,1} = mean(squeeze(dfof_align(icell, idx, range_resp)),2);
    dfof_tg_noad{icell,1} = resp_cond{icell,1} - base_cond{icell,1};
    
    mu = mean(dfof_tg_noad{icell,1});
    sigma = std(dfof_tg_noad{icell,1});
    z_score_{icell,1} = (dfof_tg_noad{icell,1} - mu) / sigma;
    z_score(icell, :) = z_score_{icell};    
end

if save_flag; save dfof_tg_noad_z_score.mat z_score dfof_tg_noad; end 

end