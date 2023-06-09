function dfof_align = dfof_by_trial_base(tc_aligned, npSub_tc, frame_ad)
global frame_rate ntrial ncell

% input: 
%     tc_aligned by ad or tg. ncell x ntrial x trial_len
%     npSub_tc. nframe x ncell
%     frame_ad as end point of trial-specific baseline
%     frame_rate, ncell, ntrial
% output:
%     dfof_align = (tc - base) / base. ncell x ntrial x trial_len
% 
% TODO: fix baseline bug!

% %% rewrite 2022-04-18
% 
% trial_base_len = frame_rate * 1; % 30 frame/sec * 1 sec
% % trial_base_len = floor(frame_rate * 0.2);
% tc_trial_base = tc_aligned(:, :, end-trial_base_len:end);
% 
% base_max = max(tc_trial_base, [], 3);
% base_min = min(tc_trial_base, [], 3);
% base_fluc = (base_max - base_min) ./ (base_min + 1e-6);
% base_fluc = mean(base_fluc(:))
% if base_fluc >= 0.1; disp('trial baseline fluctuates too much!'); end
% 
% tc_trial_base = mean(tc_trial_base, 3);
% dfof_align = tc_aligned ./ tc_trial_base - 1;

%%
% trial_base_len = frame_rate * 1; % 30 frame/sec * 1 sec
% trial_base_len = frame_rate * 0.5;
trial_base_len = frame_rate * 0.1; % reduced trial base length for cellpose bunnytop
tc_trial_base = zeros(ncell, ntrial, trial_base_len);
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        tc_trial_base(icell, itrial, :) = ...
            [npSub_tc_cell(frame_ad(itrial) - trial_base_len+1 : frame_ad(itrial))];
    end
end

t = squeeze(nanmean(squeeze(tc_trial_base(:,:,:)), 1)); 
t_base = squeeze(nanmean(t(:,:), 1)); 
alarm = (max(t_base) - min(t_base)) / min(t_base);
if alarm >= 0.05; disp('trial baseline fluctuates too much!'); end
% figure;
% plot(t_base, 'k'); % saveas(gcf, ['trial base'], 'jpg'); close 

tc_trial_base_avg = nanmean(tc_trial_base, 3);
% figure
% imagesc(tc_trial_base_avg) % no overall drift of baseline across trials in session

dfof_align = zeros(ncell, ntrial, size(tc_aligned, 3));
for icell = 1:ncell
for itrial = 1:ntrial
    dfof_align(icell, itrial, :) = tc_aligned(icell, itrial, :) ./ tc_trial_base_avg(icell, itrial, :) - 1;
end
end
% 1+1