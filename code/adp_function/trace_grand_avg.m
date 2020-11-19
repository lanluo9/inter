function [trace_avg, trace_sem]= trace_grand_avg(dfof_align_ad, save_flag)

% input: dfof_align by ad or tg, linked w/ resp_win 'ad' or 'tg'. save_flag to save mat
% output: trace avg and sem aligned by ad = ncell x nori x nisi3 [noad 750 250]

global ncell nori id_ori id_isi3 trial_len_max

trace_by_trial = cell(ncell, nori, length(id_isi3)); 
trace_avg = cell(ncell, nori, length(id_isi3)); 
trace_sem = cell(ncell, nori, length(id_isi3)); 

for icell = 1 : ncell    
for iori = 1 : nori 
    for iisi =  1 : length(id_isi3) 
        idx = intersect(id_ori{iori}, id_isi3{iisi}); 
        trace_by_trial{icell, iori, iisi} = squeeze(dfof_align_ad(icell, idx, 1:trial_len_max));
        trace_avg{icell, iori, iisi} = nanmean(trace_by_trial{icell, iori, iisi}, 1);
        trace_sem{icell, iori, iisi} = nanstd(trace_by_trial{icell, iori, iisi}, 1) ./ length(idx);
    end
end
end

if save_flag; save trace_aligned.mat trace_avg trace_sem; end