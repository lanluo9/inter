function dfof_align = dfof_by_trial_base(tc_aligned, npSub_tc, frame_ad, paradigm_ms)
global frame_rate ntrial ncell

% input: 
%     tc_aligned by ad or tg. ncell x ntrial x trial_len
%     npSub_tc. nframe x ncell
%     frame_ad as end point of trial-specific baseline
%     frame_rate, ncell, ntrial
% output:
%     dfof_align = (tc - base) / base. ncell x ntrial x trial_len

% %% rewrite 2023-06-12 
% % this seems to fix baseline bug. now baseline drops back to around 0 at the end of ITI
% % but only bc the baseline here is taken from end of ITI
% 
% % half_off_period = iti_ms - 2 * stim_ms - max_isi_ms;
% % tc_aligned_shift = tc_aligned; % shift trace to align not with stim1 onset, but 
% 
% trial_base_len = frame_rate * 1; % 30 frame/sec * 1 sec
% tc_trial_base = tc_aligned(:, :, end-trial_base_len:end);
% 
% % base_max = max(tc_trial_base, [], 3);
% % base_min = min(tc_trial_base, [], 3);
% % base_fluc = (base_max - base_min) ./ (base_min + 1e-6); % for cell, for trial, get baseline fluctuation
% % base_fluc = mean(base_fluc(:))
% % if base_fluc >= 0.1; disp('trial baseline fluctuates too much!'); end
% 
% tc_trial_base = mean(tc_trial_base, 3);
% dfof_align = tc_aligned ./ tc_trial_base - 1;

%%
trial_base_len = frame_rate * 1; % 30 frame/sec * 1 sec
% trial_base_len = frame_rate * 0.1; % reduced trial base length for cellpose bunnytop
tc_trial_base = zeros(ncell, ntrial, trial_base_len);
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        tc_trial_base(icell, itrial, :) = ...
            [npSub_tc_cell(frame_ad(itrial) - trial_base_len+1 : frame_ad(itrial))];
    end
end

t = squeeze(nanmean(squeeze(tc_trial_base(:,:,:)), 1)); % agg over cells
t_base = squeeze(nanmean(t(:,:), 1)); % agg over trials. len = nframe_base
alarm = (max(t_base) - min(t_base)) / min(t_base);
if alarm >= 0.05; disp('trial baseline fluctuates too much!'); end
% figure;
% plot(t_base, 'k'); % saveas(gcf, ['trial base'], 'jpg'); close 

tc_trial_base_avg = nanmean(tc_trial_base, 3); % agg over nframe_base
% figure
% imagesc(tc_trial_base_avg) % no overall drift of baseline across trials in session

dfof_align = zeros(ncell, ntrial, size(tc_aligned, 3));
for icell = 1:ncell
for itrial = 1:ntrial
    dfof_align(icell, itrial, :) = tc_aligned(icell, itrial, :) ./ tc_trial_base_avg(icell, itrial, :) - 1;
end
end

%% plot trial trace as: off/2 -> on -> off/2 
% %% dfof_by_trial_base check
% plot from previous trial's stim2 offset to next trial's stim1 onset
% a full off-on-off cycle, to see if baseline falls back properly

% iti_ms = 6000; % TODO: dont hard code it
% stim_ms = 100;
% max_isi_frame = 23;
% max_isi_ms = max_isi_frame / frame_rate * 1000;
half_off_ms = (paradigm_ms.iti_ms...
               - paradigm_ms.stim1_ms - paradigm_ms.stim2_ms ...
               - paradigm_ms.max_isi_ms) / 2;
half_off_frame = frame_rate * half_off_ms/1000;
half_off_frame = int64(half_off_frame);

dfof_align_shift = dfof_align(:, :, 1:(end-half_off_frame)); % on -> off/2
dfof_align_tail = dfof_align(:, :, (end-half_off_frame):end); % off/2

dfof_align_tail(:, end, :) = 0; % final ITI/2 is useless
tmp = dfof_align_tail(:, 1:(end-1), :);
dfof_align_tail = cat(2, dfof_align_tail(:, end, :), tmp); % move final ITI/2 to first trial to act as padding

dfof_align_shift = cat(3, dfof_align_tail, dfof_align_shift);
% size(dfof_align_shift) % ncell x ntrial x nframe

t1 = squeeze(nanmean(dfof_align_shift, 1));
t1 = squeeze(nanmean(t1, 1)); 
t2 = squeeze(nanmedian(dfof_align_shift, 1));
t2 = squeeze(nanmedian(t2, 1)); 

figure
subplot(121)
plot(t1)
hold on
yline(0)
yline(0.01)
xline(double(half_off_frame))
title('mean')

subplot(122)
plot(t2)
hold on
yline(0)
yline(0.01)
xline(double(half_off_frame))
title('median')

