function [ppResp] = dfof_resp_trialwise_jeff(dfof_align_tg, save_flag)
% write data in correct format for jeff population vector decoder - jin2019
% input: 
    % always use dfof_align_tg for this use case, bc we only need R2 (adapted response)
% output:
    % ppResp = "paired pulse response"
    % [nisi x nori] matlab-cells, where 
    % isi order = 250-750-inf
    % ori order = 22-90-157-0 deg = stim2 ori id (1-7-0)
    % in each matlab-cell = [ncell x ntrial] mat

global ncell nori nisi id_isi3 id_ori range_base range_resp
id_isi3_jeff = id_isi3([3, 2, 1]); % original isi order is inf-750-250) need to reverse to 250-750-inf
id_isi3_jeff = id_isi3_jeff([1, 3, 2]); % want to compare isi 6k to 250, so re-order isi to 250-inf-750
id_ori_jeff = id_ori([2; 3; 4; 5; 6; 7; 8; 1]); % change ori order to: 0 deg (180 deg) at the end

ppResp = cell(nisi, nori); 

for iisi =  1 : length(id_isi3_jeff) 
for iori = 1 : nori 
    trial_id = intersect(id_ori_jeff{iori}, id_isi3_jeff{iisi}); 

    base_cond = squeeze(mean(dfof_align_tg(:, trial_id, range_base), 3)); 
    resp_cond = squeeze(mean(dfof_align_tg(:, trial_id, range_resp), 3));
    ppResp{iisi, iori} = resp_cond - base_cond;

end
end

if save_flag; save dfof_trial_jeff.mat ppResp; end

end