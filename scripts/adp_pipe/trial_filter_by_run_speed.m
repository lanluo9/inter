%% path

mouse = 'i1375';
date = '220915';
run = '003';
time = '1421';
area = 'V1';

%% load mworks data

fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
input_mworks = load(fName);

speed = input_mworks.input.wheelSpeedValues;
assert(length(speed) == input_mworks.input.trialSinceReset)
size(speed)
size(speed{1,1})
size(speed{1,2})

%% align speed

trial_len_arr = cellfun(@length, speed);
trial_len_min = min(trial_len_arr)

speed_trim = cellfun(@(v) v(1:trial_len_min), speed, 'uniformoutput', false); % trim trials to min len
speed_trim = cell2mat(cellfun(@transpose, speed_trim, 'uniform', 0));
speed_trim = speed_trim'; % ntrial x nframe
ntrial = size(speed_trim, 1)

%% visliz

figure;
imagesc(speed_trim); % some trials have higher running speed. no frames in trial have very high run speed
colorbar;

speed_avg = mean(speed_trim, 1); % run speed of avg trial
figure;
plot(speed_avg);

speed_dist = mean(speed_trim, 2); % speed dist across trials
figure;
histogram(speed_dist, 500);

%% discard trial by run speed

discard_perc_low = 1
discard_perc_high = 4.75
thres_low = prctile(speed_dist, discard_perc_low)
thres_high = prctile(speed_dist, 100-discard_perc_high)

figure;
histogram(speed_dist, 200);
hold on;
xline(thres_low);
xline(thres_high);

ntrial_discard = sum((speed_dist < thres_low) | (speed_dist > thres_high))
ntrial_discard / ntrial % percent of discarded trials should be close to discard_perc_low + high
ntrial - ntrial_discard

trial_speed_discard = (speed_dist < thres_low) | (speed_dist > thres_high); % discard trials with extreme speed
trial_speed_ok = ~trial_speed_discard;

%% visliz after discard

figure;
subplot(1,2,1)
imagesc(speed_trim);
subplot(1,2,2)
imagesc(speed_trim(~trial_speed_discard, :))
colorbar;
set(gcf, 'Position', get(0, 'Screensize'));

figure;
speed_avg = mean(speed_trim, 1);
plot(speed_avg);
hold on;
speed_avg_filtered = mean(speed_trim(~trial_speed_discard, :), 1);
plot(speed_avg_filtered);

%% save trial filter with cell filter

mat_inter_path = ['Z:\All_Staff\home\lan\Data\2P_images\mat_inter\', area, '_', mouse, '_', date, '_cellpose']
cd(mat_inter_path)
speed_trial_by_frame = speed_trim;
save('filter_trials_by_speed.mat', 'trial_speed_ok', 'speed_trial_by_frame')




