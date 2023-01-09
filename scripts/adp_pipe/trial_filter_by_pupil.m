%%

close all
clear all
clc

mouse = 'i1375';
date = '220915';
run = '003'; % mix50 session 1
time = '1421';
area = 'V1';
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
CD = [data_pn '\Data\2P_images\' mouse '\' date '\' run]; % Z:\All_Staff\home\lan\Data\2P_images\i1375\220915\003
cd(CD);
fn = [run '_000_000_eye.mat'];

%%

% load data
data_temp = load(fn);
data_temp = squeeze(data_temp.data);

% crop frames to match mworks data
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nFrames = input.counterValues{end}(end);
data = data_temp(:,:,1:nFrames);      % the raw images...

% data = data(:, :, 1:1000); 
% disp('using test data only')

%% Crop image to isolate pupil 

% select an area as small as possible that contains pupil, but excludes bright spots
[data_crop, rect] = cropEyeData(data); % select rectangle, hit enter if happy, z to redo

%% measure pupil position/diameter

rad_range = [3 13]; % histogram of pupil radius must show both tails, otherwise adjust rad_range
Eye_data = extractEyeData(data_crop, rad_range); % check if pupil not found in too many frames
% must use full data (not truncated test data) below, to match mworks input

%% align to stimulus presentation

[rad, centroid] = alignEyeData(Eye_data, input); % eye data peri-stimulus 1
[rad2, centroid2] = alignEyeData(Eye_data, input, 2); % eye data peri-stimulus 2

%% visliz radius

rad_avg = (mean(rad.tc, 2) + mean(rad2.tc, 2)) / 2;
figure
plot(rad_avg)
hold on
xline(3.5)
xlabel('pupil radius avg over trials')

figure
scatter(rad.base, rad.stim, 'MarkerEdgeAlpha',.4);
hold on
scatter(rad2.base, rad2.stim, 'MarkerEdgeAlpha',.4);
plot(xlim,ylim,'-k')
xlabel('baseline pupil radius')
ylabel('stimulus on pupil radius') % stim-on pupil slightly larger than baseline

rad_trials = (mean(rad.tc, 1) + mean(rad2.tc, 1)) / 2;
figure;
histogram(rad_trials); % should discard trials w extreme pupil size

%% visliz position

% figure
% scatter(centroid.base(1, :), centroid.base(2, :));
% hold on
% scatter(centroid.stim(1, :), centroid.stim(2, :));

% ntrial = length(centroid.dist);
% figure
% hold on
% for itrial = 1 : ntrial
%     plot( [centroid.base(1, itrial), centroid.stim(1, itrial)], ...
%         [centroid.base(2, itrial), centroid.stim(2, itrial)] )
% end % should discard trials w large diff between baseline and stim-on
% xlabel('eye position, moving from baseline to stim-on in each trial')

% figure % same as lineplot above but with arrow
% hold on
% for itrial = 1 : ntrial
%     p1 = [centroid.base(1, itrial), centroid.base(2, itrial)];
%     p2 = [centroid.stim(1, itrial), centroid.stim(2, itrial)];
%     dp = p2 - p1; % position difference
%     quiver(p1(1),p1(2), dp(1),dp(2), 0)
% end

ntrial = length(centroid.dist);
figure
hold on
for itrial = 1 : ntrial % same as above, but for eye position distance btw stim 1 vs 2 
    plot( [centroid.stim(1, itrial), centroid2.stim(1, itrial)], ...
        [centroid.stim(2, itrial), centroid2.stim(2, itrial)] )
end % should discard trials w large diff 
xlabel('eye position, moving from stim1 to stim2 in each trial')

dp = zeros(ntrial, 2);
for itrial = 1 : ntrial
    p1 = [centroid.stim(1, itrial), centroid.stim(2, itrial)];
    p2 = [centroid2.stim(1, itrial), centroid2.stim(2, itrial)];
    dp(itrial, :) = p2 - p1; % position difference
end
eye_move_dist = sqrt(dp(:, 1).^2 + dp(:, 2).^2);
figure;
histogram(eye_move_dist);
xlabel('eye movement distance from stim1 to stim2 within each trial')

pupil_deviation = (centroid.dist + centroid2.dist) / 2;
figure;
histogram(pupil_deviation); % should discard trials w pupils far from median position
xlabel('pupil centroid distance from median')

%% discard trial by pupil radius & position

discard_perc_low = 1
discard_perc_high = 5

% extreme pupil size
[rad_thres, rad_id_retained]= threshold_percentile(rad_trials, discard_perc_low, discard_perc_high);

% large diff between stim 1 vs stim 2
[eyemove_thres, eyemove_id_retained]= threshold_percentile(eye_move_dist, 0, discard_perc_high);

% pupils far from median position
[pupil_deviation_thres, pupil_deviation_id_retained]= threshold_percentile(pupil_deviation, 0, discard_perc_high);

%% visliz after discard

figure;
histogram(rad_trials);
hold on
histogram(rad_thres);

figure;
histogram(eye_move_dist);
hold on
histogram(eyemove_thres);

figure;
histogram(pupil_deviation);
hold on
histogram(pupil_deviation_thres);

trial_eye_ok = (rad_id_retained & eyemove_id_retained & pupil_deviation_id_retained);
sum(trial_eye_ok) / length(trial_eye_ok) % ratio of trials that passed pupil check

mat_inter_path = ['Z:\All_Staff\home\lan\Data\2P_images\mat_inter\', area, '_', mouse, '_', date, '_cellpose']
cd(mat_inter_path)
save('filter_trials_by_pupil.mat', 'trial_eye_ok', 'rad_trials', 'eye_move_dist', 'pupil_deviation')


