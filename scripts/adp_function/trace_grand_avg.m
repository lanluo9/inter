function [trace_avg, trace_sem]= trace_grand_avg(dfof_align_ad, save_flag)

% input: dfof_align by ad or tg, linked w/ resp_win 'ad' or 'tg'. save_flag to save mat
% output: trace avg and sem aligned by ad = ncell x nori x nisi3 [noad 750 250]
% 2021-03-29 update: trace = ncell x nori x nisi3 x nframe_trial (padded with NaN)
% 2021-04-02 update: use min trial len to cut off the tail, instead of padding NaN

global ncell nori id_ori id_isi3 trial_len_min

trace_avg = zeros(ncell, nori, length(id_isi3), trial_len_min); 
trace_sem = zeros(ncell, nori, length(id_isi3), trial_len_min); 

% for icell = 1 : ncell    
for iori = 1 : nori 
    % disp(['ori # ', num2str(iori), ' of 8'])
    for iisi =  1 : length(id_isi3) 
        % disp(['  isi # ', num2str(iisi), ' of 3'])
        idx = intersect(id_ori{iori}, id_isi3{iisi}); 
        temp = dfof_align_ad(:, idx, 1:trial_len_min); % should not squeeze, leave idx (axis=2) to be averaged
        trace_avg(:, iori, iisi, :) = squeeze(nanmean(temp, 2));
        trace_sem(:, iori, iisi, :) = squeeze(nanstd(temp, 2)) ./ length(idx);
    end
end
% end

% trace_avg = squeeze(trace_avg);

if save_flag; save trace_aligned.mat trace_avg trace_sem; end
