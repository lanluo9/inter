%%

close all
clear all
clc

mouse = 'i1372';
date = '220714';
area = 'V1';

data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
data_path = [data_path '\Data\2P_images\' mouse '\' date]; 
cd(data_path);

master_xls = '2P Imaging Notes Lan.xlsx';
dataset_meta = readtable(master_xls);
dataset_meta = dataset_meta(dataset_meta.Var7 > 10000, :); % nframes large enough to be grat_6SF session

nset = size(dataset_meta,1);

%%
for i = 1 : nset
    run = dataset_meta{i,1}{1}(1:3)
    time = num2str(dataset_meta{i,8})

%% load pupil data
cd([data_path, '\', run]);
fn = [run '_000_000_eye.mat'];
data_temp = load(fn);
data_temp = squeeze(data_temp.data);

% crop frames to match mworks data
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nFrames = input.counterValues{end}(end);
data = data_temp(:,:,1:nFrames);      % the raw images...

data_test = data(:, :, 1:1000); 
disp('using test data only')

%% Crop image to isolate pupil 

% select an area as small as possible that contains pupil, but excludes bright spots
[data_crop, rect] = cropEyeData(data_test); % select rectangle, hit enter if happy, z to redo

%% measure pupil position/diameter

rad_range = [2 13]; % histogram of pupil radius must show both tails, otherwise adjust rad_range
Eye_data = extractEyeData(data_crop, rad_range); % check if pupil not found in too many frames
% must use full data (not truncated test data) below, to match mworks input

%% rerun w full eye data

% pause % TODO: add user input? to adjust rad range after test data
close all
data_crop = data(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:); % crop full data by rectangle determined above
Eye_data = extractEyeData(data_crop, rad_range);

%% align to stimulus presentation

close all
[rad, centroid] = alignEyeData(Eye_data, input); % eye data peri-stimulus 1
[rad2, centroid2] = alignEyeData(Eye_data, input, 2); % eye data peri-stimulus 2

%% visliz radius

close all
rad_avg = (nanmean(rad.tc, 2) + nanmean(rad2.tc, 2)) / 2;
figure
plot(rad_avg)
hold on
xline(3.5)
xlabel('pupil radius avg across frames of a trial')

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
xlabel('pupil size distribution across trials')

%% visliz position

close all

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

close all;

figure;
histogram(rad_trials);
hold on
histogram(rad_thres);
xlabel('pupil size distribution across trials')

figure;
histogram(eye_move_dist);
hold on
histogram(eyemove_thres);
xlabel('eye movement distance stim1-2. distribution across trials')

figure;
histogram(pupil_deviation);
hold on
histogram(pupil_deviation_thres);
xlabel('pupil position deviation from median distribution across trials')

trial_eye_ok = (rad_id_retained & eyemove_id_retained & pupil_deviation_id_retained);
disp('ratio of trials that passed pupil check')
sum(trial_eye_ok) / length(trial_eye_ok) % ratio of trials that passed pupil check

%% load speed data

close all
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
input_mworks = load(fName);

speed = input_mworks.input.wheelSpeedValues;
assert(length(speed) == input_mworks.input.trialSinceReset)
% size(speed)
% size(speed{1,1})
% size(speed{1,2})

%% align speed

trial_len_arr = cellfun(@length, speed);
trial_len_min = min(trial_len_arr)

speed_trim = cellfun(@(v) v(1:trial_len_min), speed, 'uniformoutput', false); % trim trials to min len
speed_trim = cell2mat(cellfun(@transpose, speed_trim, 'uniform', 0));
speed_trim = speed_trim'; % ntrial x nframe
speed_trial_by_frame = speed_trim;
ntrial = size(speed_trim, 1);

%% visliz

speed_flat = reshape(speed_trim, [size(speed_trim,1)*size(speed_trim,2),1]);

figure;
imagesc(speed_trim); % some trials have higher running speed. no frames in trial have very high run speed
colorbar;

speed_avg = mean(speed_trim, 1); % run speed of avg trial
figure;
plot(speed_avg);

speed_dist = mean(speed_trim, 2); % avg speed dist across trials
figure;
histogram(speed_dist, 500);

%% discard trial by run speed

close all

discard_perc_low = 1
discard_perc_high = 4.75
[speed_thresholded, trial_speed_ok] = threshold_percentile(speed_dist, discard_perc_low, discard_perc_high);

% figure;
% histogram(speed_dist, 200);
% hold on;
% xline(thres_low);
% xline(thres_high);

%% visliz after discard

close all

figure;
subplot(1,2,1)
imagesc(speed_trim);
subplot(1,2,2)
imagesc(speed_trim(trial_speed_ok, :))
colorbar;
set(gcf, 'Position', get(0, 'Screensize'));

figure;
speed_avg = mean(speed_trim, 1);
plot(speed_avg);
hold on;
speed_avg_filtered = mean(speed_trim(trial_speed_ok, :), 1);
plot(speed_avg_filtered);

%% save trial filter with cell filter

close all

mat_inter_path = ['Z:\All_Staff\home\lan\Data\2P_images\mat_inter\', area, '_', mouse, '_', date, '_cellpose']
try cd(mat_inter_path)
catch
    mkdir(mat_inter_path)
    cd(mat_inter_path)
end

save(['filter_trials_by_pupil_or_speed_', run, '.mat'], 'trial_speed_ok', 'speed_trial_by_frame', ...
    'trial_eye_ok', 'rad_trials', 'eye_move_dist', 'pupil_deviation')

end