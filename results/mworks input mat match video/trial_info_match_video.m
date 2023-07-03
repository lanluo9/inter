%%

dir_mat = 'C:\Users\ll357\Documents\inter\results\mworks input mat match video\data-i9999-230630-1121.mat'
tmp = load(dir_mat);
input_behav = tmp.input;
clear tmp

contrast_ad = celleqel2mat_padded(input_behav.tBaseGratingContrast)'; 
% sum(contrast_ad) / length(contrast_ad)

frame_ad = double(cell2mat(input_behav.cStimOn)); frame_ad_off = double(cell2mat(input_behav.cStimOff));
frame_tg = celleqel2mat_padded(input_behav.cTargetOn); frame_tg = double(frame_tg);
isi_seq = (frame_tg - frame_ad_off)';
% unique(isi_seq)

ori_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg)'; 
ori_seq(ori_seq == 180) = 0;
ori_seq = uint8(ori_seq);
% unique(ori_seq)

trial_info = [contrast_ad, isi_seq, ori_seq];

% frame_rate = 30;
% paradigm_ms.stim1_ms = input_behav.stimOnTimeMs;
% paradigm_ms.stim2_ms = input_behav.targetOnTimeMs;
% paradigm_ms.max_isi_ms = max(isi_seq) / frame_rate * 1000;
% paradigm_ms.iti_ms = input_behav.itiTimeMs;








%%

% tmp = squeeze(nanmean(nanmean(dfof_align_ad, 1), 2));
% % tmp = squeeze(dfof_align_ad(25, 16, :));
% figure
% plot(tmp)
% hold on
% yline(0)
% % xlim([0, 60]);
% 
% % 10(11), 20(19), 34
% % 8 or 10 frame = 250 ms
% % 24 frame = 750 ms
% 
% %%
% 
% figure
% imagesc(tc_trial_base_avg)
% 
% %%
% 
% dfof_align_ad_respwin = squeeze(nanmean(dfof_align_ad(:, :, 20:30), 3));
% figure; imagesc(dfof_align_ad_respwin)
% colorbar