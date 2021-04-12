function tc_align = align_tc(landmark, npSub_tc)
global ntrial ncell

% input: 
%     "landmark" = align by adapter or target. array of frame numbers
%     "npSub_tc" = timecourse
%      ncell, ntrial
% 
% output:
%     tc_align = aligned by adapter or target. ncell x ntrial x trial_len

trial_len = diff(landmark);
tc_align = zeros(ncell, ntrial, min(trial_len));

% for icell = 1:ncell
%     npSub_tc_cell = npSub_tc(:,icell);
%     
%     for itrial = 1:ntrial
%         start_id = landmark(itrial);
%         tc_align(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + min(trial_len) - 1)];
% %         tc_align(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + trial_len(itrial) - 1); ...
% %             NaN(max(trial_len) - trial_len(itrial), 1)];
%     end
%     
% end

for itrial = 1:ntrial
%     if mod(itrial,10) == 0; disp([num2str(itrial),' trial out of ', num2str(ntrial)]); end
    start_id = landmark(itrial);
    tc_align(:, itrial, :) = npSub_tc(start_id : (start_id + min(trial_len) - 1), :)';
end