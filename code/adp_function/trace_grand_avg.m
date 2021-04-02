function [trace_avg, trace_sem]= trace_grand_avg(dfof_align_ad, save_flag)

% input: dfof_align by ad or tg, linked w/ resp_win 'ad' or 'tg'. save_flag to save mat
% output: trace avg and sem aligned by ad = ncell x nori x nisi3 [noad 750 250]
% 2021-03-29 update: trace = ncell x nori x nisi3 x nframe_trial (padded with NaN)

global ncell nori id_ori id_isi3 trial_len_max

trace_avg = zeros(ncell, nori, length(id_isi3), trial_len_max); 
trace_sem = zeros(ncell, nori, length(id_isi3), trial_len_max); 

% for icell = 1 : ncell    
for iori = 1 : nori 
    disp(iori)
    for iisi =  1 : length(id_isi3) 
        disp(iisi)
        idx = intersect(id_ori{iori}, id_isi3{iisi}); 
        temp = squeeze(dfof_align_ad(:, idx, 1:trial_len_max));
        trace_avg(:, iori, iisi, :) = nanmean(temp, 1);
        trace_sem(:, iori, iisi, :) = nanstd(temp, 1) ./ length(idx);
    end
end
% end

if save_flag; save trace_aligned.mat trace_avg trace_sem; end
