function [data_reg, LL_base, date, imouse, run_str] = get_data_reg_cellpose_tif(...
    arg_mouse, arg_date, arg_ImgFolder, stim_type, run_str_ref)

mouse = arg_mouse;
[data, behav_input, LL_base, date, imouse, run_str] = load_sbx_data( ...
    arg_mouse, arg_date, arg_ImgFolder);

% check data corruption
data_frame_avg = squeeze(mean(mean(data, 1) ,2));
frame_corrupt = min(find(data_frame_avg == 65535)) - 1; % fully corrupted frame avg = 65535
if ~isempty(find(data_frame_avg == 65536))
    sprintf('data corruption starts at frame no. %d', frame_corrupt)

    nep = floor(size(data,3)./10000);
    [n, n2] = subplotn(nep);
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nep
        subplot(n,n2,i); 
        imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); 
        title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); 
    end
    print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_check_corrupt.pdf']),'-dpdf', '-bestfit')
else
    sprintf('data corruption not detected')
end

% Choose register interval
if strcmp(run_str_ref, run_str) % if this session is ref session, then register against oneself
    nep = floor(size(data,3)./10000);
    [n, n2] = subplotn(nep);
    figure('units','normalized','outerposition',[0 0 1 1]);
    fluo_arr = zeros(nep, 1);
    for i = 1:nep
        subplot(n,n2,i); 
        imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); 
        frame_avg = mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3);
        fluo_arr(i) = mean(frame_avg(:));
        title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); 
    end
    print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_select_register_ref.pdf']),'-dpdf', '-bestfit')
%     close all
    
    disp('start registration. using middle section as motion correct ref')
    fluo_arr = fluo_arr(2:end-1); % dont use beginning or end of recording as data_avg
    [~, fluo_median_id] = min(abs(fluo_arr - median(fluo_arr))); % find arr element closest to median fluo
    select = fluo_median_id + 1 % add back the beginning
    start_idx = select * 10000 + 1;
    stop_idx = select * 10000 + 500;
    data_avg = mean(data(:,:,start_idx:stop_idx),3); % use the session itself as data_avg
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str_ref], [date '_' imouse '_' run_str_ref '_data_avg.mat']), 'data_avg')
else % otherwise take data_avg from ref session
    try
        load(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str_ref], [date '_' imouse '_' run_str_ref '_reg_shifts.mat']), 'data_avg')
    catch
        load(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str_ref], [date '_' imouse '_' run_str_ref '_data_avg.mat']), 'data_avg')
    end
end

close all
% Register data
if exist(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']))
    disp('re-register from reg_shifts.mat')
    load(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'behav_input')
    [outs, data_reg]=stackRegister_MA_LL(double(data),[],[],out);
% for grating dataset 2 - i1323 V1:
% Index in position 1 exceeds array bounds (must not exceed 99865).
% 
% Error in stackRegister_MA (line 76)
%         row_shift = shifts_xy(index,3);
% 
% Error in get_data_reg_cellpose_tif (line 73)
%     [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
% 
% Error in batch_cellpose (line 55)
% [data_reg, LL_base, date, imouse, run_str] =
% get_data_reg_cellpose_tif(...
    clear out outs
else
    tic; [out, data_reg] = stackRegister_LL(data, data_avg); toc;
%     mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'behav_input')
end
clear data out
disp('registration finished')

% test stability
nep = floor(size(data_reg,3)./10000);
[n, n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%%

disp(stim_type)
disp(run_str)

switch stim_type % TODO: wrap into function
case {'bunny', 'mix', 'grat6'}
    
%% find activated cells

cStart = cell2mat(behav_input.cStimOneOn); % same as cStimOn
cStimOn = cell2mat(behav_input.cStimOneOn);
cStimOff = cell2mat(behav_input.cStimOneOff);
cTarget = cell2mat(behav_input.cStimTwoOn);
try nTrials = behav_input.trialsSinceReset % 464 = 32 types * 14.5 reps
catch
    nTrials = behav_input.trialSinceReset
end
sz = size(data_reg); % [y pixel * x pixel * nframe]
data_f = zeros(sz(1),sz(2),nTrials);
data_adapter = zeros(sz(1),sz(2),nTrials);
data_f2 = zeros(sz(1),sz(2),nTrials);
data_targ = zeros(sz(1),sz(2),nTrials);

% show trial design
unique(cStimOff - cStimOn) % adapter = 3-4 frame. rounded from 100 ms?
unique(cTarget - cStimOff) % ISI = 8 or 22-23 frame. rounded from 250? or 750 ms
unique(cTarget - cStimOn) % target: 3-4 frame too
trial_len = diff(cStart); unique(trial_len)

% determine ca signal latency
tc_screen = mean(mean(data_reg,1),2);
tc_screen = squeeze(tc_screen);
whos tc_screen
data_trial = zeros(min(unique(trial_len)), nTrials); % take 1-200 frame of every trial
data_trial_real = zeros(max(trial_len), nTrials);
for it = 1:(nTrials-1)
    start_id = cStimOn(it);
    data_trial(:,it) = tc_screen(start_id : start_id + min(unique(trial_len)) - 1);
    data_trial_real(:,it) = [tc_screen(start_id : start_id + trial_len(it) - 1); NaN(max(trial_len) - trial_len(it), 1)];
end

figure; plot(mean(data_trial, 2))
set(gcf, 'Position', get(0, 'Screensize'));
analysis_dir = fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]);
cd(analysis_dir)
saveas(gcf, ['find_ca_latency.jpg'])
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_find_ca_latency.pdf']),'-dpdf', '-bestfit')
close

figure; data_trial_zoom_in = nanmean(data_trial_real, 2); plot(data_trial_zoom_in(1:40)); grid on; grid minor
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['find_ca_latency_zoomin.jpg'])
save find_ca_latency.mat data_trial_zoom_in data_trial
disp('detecting first peak automatically')
disp(' ')

%% calculate response

find_peak_bef = 16;
disp('assume: first peak comes before 16 frames. second peak comes after')
data_trial_zoom_first_peak = data_trial_zoom_in(1:find_peak_bef);
[~, peak_id] = max(data_trial_zoom_first_peak)
if peak_id < 6
    disp('WARNING: strange trace or peak')
end

window_len = 2; % 2-3
ca_latency = peak_id - floor(window_len/2) - 1; % - half window - starting frame number (1)
% ca_latency = x-1. stim onset frame 1 -> signal received frame x

figure;
plot(data_trial_zoom_in(1:40));
xline(1 + ca_latency, 'k--')
xline(1 + ca_latency + window_len, 'k--')
grid on; grid minor
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['find_ca_latency_ca_window.jpg'])
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_find_ca_latency_ca_window.pdf']),'-dpdf', '-bestfit')
close all

assert(length(cTarget) == nTrials && length(cStart) == nTrials && cTarget(nTrials)+3 < sz(3))
for itrial = 1:nTrials
        data_f(:,:,itrial) = mean(data_reg(:,:,...
            cStart(itrial) - 10 : cStart(itrial) - 1),3);
        
        data_adapter(:,:,itrial) = mean(data_reg(:,:,...
            cStimOn(itrial) + ca_latency : cStimOn(itrial) + ca_latency + window_len),3);
        
        data_f2(:,:,itrial) = mean(data_reg(:,:,...
            cTarget(itrial) - 3 :cTarget(itrial) - 1),3); % cTarget(itrial) + 14 :cTarget(itrial) + 16),3);

        data_targ(:,:,itrial) = mean(data_reg(:,:,...
            cTarget(itrial) + ca_latency : cTarget(itrial) + ca_latency + window_len),3);
end

data_adapter_dfof = (data_adapter-data_f)./data_f;
% data_base2_dfof = (data_f2-data_f)./data_f;
data_targ_dfof = (data_targ-data_f2)./data_f2; 
data_targ_dfof_fake = (data_targ-data_f)./data_f; 

%% plot dfof for randStim1_doSameStim2 bunny

% assert(behav_input.doRandStimOne == 1 & behav_input.doSameStims == 1)
% adapter_id = cell2mat(behav_input.tstimOne);
% adapter_list = unique(adapter_id);
% n_adapter = length(adapter_list);
% target_id = cell2mat(behav_input.tstimTwo);
% target_list = unique(target_id);
% n_target = length(target_list);

% try
%     assert(behav_input.doRandStimOne == 1 & behav_input.doSameStims == 1) % bunny 
% catch
%     assert(behav_input.doRandSF == 1) % or grat_SF6
% end

try % randStim1_doSameStim2 with bunnies6.mwel
    adapter_id = cell2mat(behav_input.tstimOne);
    target_id = cell2mat(behav_input.tstimTwo);
catch % grat_SF6 with twoStim.mwel
    SF_arr = sort(unique(cell2mat(behav_input.tStimOneGratingSpatialFreqCPD)));
    adapter_SF = cell2mat(behav_input.tStimOneGratingSpatialFreqCPD);
    adapter_id = zeros(length(adapter_SF), 1);
    for i = 1 : length(adapter_SF)
        adapter_id(i) = find(SF_arr==adapter_SF(i));
    end
    tmp = cell2mat(behav_input.tStimTwoGratingSpatialFreqCPD) == cell2mat(behav_input.tStimOneGratingSpatialFreqCPD);
    assert(sum(tmp) == length(tmp)) % stim 1 vs 2 have same SF
    target_id = adapter_id;
end
adapter_list = unique(adapter_id);
n_adapter = length(adapter_list);
target_list = unique(target_id);
n_target = length(target_list);

data_dfof_ad = zeros(sz(1),sz(2),n_adapter);
for i_ad = 1:n_adapter
    ind = find(adapter_id == adapter_list(i_ad));
    data_dfof_ad(:,:,i_ad) = nanmean(data_adapter_dfof(:,:,ind),3);
end

data_dfof_targ = zeros(sz(1),sz(2),n_target);
data_dfof_targ_fake = zeros(sz(1),sz(2),n_target);
for i_tg = 1:n_target
    ind = find(target_id == target_list(i_tg));
    data_dfof_targ(:,:,i_tg) = nanmean(data_targ_dfof(:,:,ind),3);
    data_dfof_targ_fake(:,:,i_tg) = nanmean(data_targ_dfof_fake(:,:,ind),3);
end

myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof = cat(3, data_dfof_ad, data_dfof_targ, data_dfof_targ_fake);
data_dfof_max = max(imfilter(data_dfof, myfilter),[],3); % gauss smooth each stim resp, then take max
data_dfof = cat(3,data_dfof, data_dfof_max); % adapter, targ, targ_fake, gaussian filter max proj

figure; imagesc(data_dfof_max)
set(gcf, 'Position', get(0, 'Screensize'));
title('data dfof max')
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_data_dfof_max.pdf']),'-dpdf', '-bestfit')
close all

    
%%
case 'grating'

% %% find activated cells

cStart = cell2mat(behav_input.cStimOn); % same as cStimOn
cStimOn = cell2mat(behav_input.cStimOn);
cStimOff = cell2mat(behav_input.cStimOff);
cTarget = cell2mat(behav_input.cTargetOn);
try nTrials = behav_input.trialsSinceReset
catch
    nTrials = behav_input.trialSinceReset
end
sz = size(data_reg); % [y pixel * x pixel * nframe]


% % for data corrupted halfway
% % data_reg = data_reg(:,:,1:87139);
% nframe_good = sz(3);
% cTarget = cTarget(cTarget <= nframe_good); % stim2 onset time must be before data corruption
% cStart = cStart(1 : length(cTarget));
% cStimOn = cStimOn(1 : length(cTarget));
% cStimOff = cStimOff(1 : length(cTarget));
% nTrials = length(cTarget);

data_f = zeros(sz(1),sz(2),nTrials);
data_adapter = zeros(sz(1),sz(2),nTrials);
data_f2 = zeros(sz(1),sz(2),nTrials);
data_targ = zeros(sz(1),sz(2),nTrials);

% show trial design
unique(cStimOff - cStimOn) % adapter = 3-4 frame. rounded from 100 ms?
unique(cTarget - cStimOff) % ISI = 8 or 22-23 frame. rounded from 250? or 750 ms
unique(cTarget - cStimOn) % target: 3-4 frame too
trial_len = diff(cStart); unique(trial_len)

% determine ca signal latency
tc_screen = mean(mean(data_reg,1),2);
tc_screen = squeeze(tc_screen);
whos tc_screen
data_trial = zeros(min(unique(trial_len)), nTrials); % take 1-200 frame of every trial
data_trial_real = zeros(max(trial_len), nTrials);
for it = 1:(nTrials-1)
    start_id = cStimOn(it);
    try
        data_trial(:,it) = tc_screen(start_id : start_id + min(unique(trial_len)) - 1);
        data_trial_real(:,it) = [tc_screen(start_id : start_id + trial_len(it) - 1); NaN(max(trial_len) - trial_len(it), 1)];
    catch
        data_trial(:,it) = 0;
        data_trial_real(:,it) = 0; % if final trial not long enough, set dfof value as 0
    end
end

figure; plot(mean(data_trial, 2))
set(gcf, 'Position', get(0, 'Screensize'));
analysis_dir = fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]);
cd(analysis_dir)
saveas(gcf, ['find_ca_latency.jpg'])
close

figure; data_trial_zoom_in = nanmean(data_trial_real, 2); 
plot(data_trial_zoom_in(1:50)); grid on; grid minor
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['find_ca_latency_zoomin.jpg'])
save find_ca_latency.mat data_trial_zoom_in data_trial
disp('detecting first peak automatically')
disp(' ')

%% calculate response

find_peak_bef = 15;
disp('assume: first peak comes before n frames. second peak comes after')
data_trial_zoom_first_peak = data_trial_zoom_in(1:find_peak_bef);
[~, peak_id] = max(data_trial_zoom_first_peak)
if peak_id < 6
    disp('WARNING: strange trace or peak')
end

window_len = 2; % 2-3
ca_latency = peak_id - floor(window_len/2) - 1; % - half window - starting frame number (1)
% ca_latency = x-1. stim onset frame 1 -> signal received frame x

figure;
plot(data_trial_zoom_in(1:40));
xline(1 + ca_latency, 'k--')
xline(1 + ca_latency + window_len, 'k--')
grid on; grid minor
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['find_ca_latency_ca_window.jpg'])
close all

% assert(length(cTarget) == nTrials && length(cStart) == nTrials && ...
%     cTarget(nTrials)+3 < sz(3))
for itrial = 1:nTrials-1 % disregard final trial, which doesnt have enough frame
        data_f(:,:,itrial) = mean(data_reg(:,:,...
            cStart(itrial) - 10 : cStart(itrial) - 1),3);
        
        data_adapter(:,:,itrial) = mean(data_reg(:,:,...
            cStimOn(itrial) + ca_latency : cStimOn(itrial) + ca_latency + window_len),3);
        
        data_f2(:,:,itrial) = mean(data_reg(:,:,...
            cTarget(itrial) - 3 :cTarget(itrial) - 1),3); % cTarget(itrial) + 14 :cTarget(itrial) + 16),3);

        data_targ(:,:,itrial) = mean(data_reg(:,:,...
            cTarget(itrial) + ca_latency : cTarget(itrial) + ca_latency + window_len),3);
end

data_adapter_dfof = (data_adapter-data_f)./data_f; % adapter - baseline 1 (before trial onset)
data_base2_dfof = (data_f2-data_f)./data_f;        % baseline 2 - baseline 1
data_targ_dfof = (data_targ-data_f2)./data_f2;     % target - baseline 2
data_targ_dfof_fake = (data_targ-data_f)./data_f;  % target - baseline 1

% %% plot response
targCon = celleqel2mat_padded(behav_input.tGratingContrast);
unique(targCon); % target contrast 1
if behav_input.doRandCon
    adapterCon = ones(size(targCon));
else
    adapterCon = celleqel2mat_padded(behav_input.tBaseGratingContrast);
end
% adapterCon = adapterCon(1:nTrials);
unique(adapterCon); % adapter contrast 0 or 1
ind_con = intersect(find(targCon == 1),find(adapterCon == 0));

adapterDir = celleqel2mat_padded(behav_input.tBaseGratingDirectionDeg);
% adapterDir = adapterDir(1:nTrials);
dirs = unique(adapterDir);
ndir = length(dirs);
targetDelta = celleqel2mat_padded(behav_input.tGratingDirectionDeg);
% targetDelta = targetDelta(1:nTrials);
deltas = unique(targetDelta);
nDelta = length(deltas);

data_dfof_dir = zeros(sz(1),sz(2),ndir);
data_dfof2_dir = zeros(sz(1),sz(2),ndir);
[n, n2] = subplotn(ndir);
% figure;
for idir = 1:ndir
    ind = setdiff(find(adapterDir == dirs(idir)),ind_con);
    data_dfof_dir(:,:,idir) = nanmean(data_adapter_dfof(:,:,ind),3);
    data_dfof2_dir(:,:,idir) = nanmean(data_base2_dfof(:,:,ind),3);
%     subplot(n,n2,idir)
%     imagesc(data_dfof_dir(:,:,idir))
%     title(dirs(idir))
end

figure
imagesc(data_dfof_dir) % adapter (deg == 0) resp is ok
title('data dfof dir')
set(gcf, 'Position', get(0, 'Screensize'));
figure
imagesc(data_dfof2_dir) % but baseline2 shows bright cells, almost same as adapter resp; bc fluores decay is slow
title('data dfof2 dir')
set(gcf, 'Position', get(0, 'Screensize'));

if sum(~isnan(data_dfof2_dir))>1
    data_dfof_dir_all = cat(3, data_dfof_dir, data_dfof2_dir);
else
    data_dfof_dir_all = data_dfof_dir;
end

data_dfof_targ = zeros(sz(1),sz(2),nDelta);
data_dfof_targ_fake = zeros(sz(1),sz(2),nDelta);
[n, n2] = subplotn(nDelta);
figure;
for idir = 1:nDelta
    ind = find(targetDelta == deltas(idir));
    data_dfof_targ(:,:,idir) = nanmean(data_targ_dfof(:,:,ind),3);
    data_dfof_targ_fake(:,:,idir) = nanmean(data_targ_dfof_fake(:,:,ind),3);

    subplot(n,n2,idir)
%     imagesc(data_dfof_targ(:,:,idir)) % targ resp shows dim and blurry cells
    imagesc(data_dfof_targ_fake(:,:,idir))
    title(deltas(idir))
end
set(gcf, 'Position', get(0, 'Screensize'));

data_dfof_targ_noadapt = zeros(sz(1),sz(2),nDelta); % TODO: why no-ad targ resp blurry?
[n, n2] = subplotn(nDelta);
figure;
for idir = 1:nDelta
    ind = find(targetDelta == deltas(idir) & adapterCon == 0);
    data_dfof_targ_noadapt(:,:,idir) = nanmean(data_targ_dfof(:,:,ind),3);

    subplot(n,n2,idir)
    imagesc(data_dfof_targ_noadapt(:,:,idir))
    title(deltas(idir))
end
set(gcf, 'Position', get(0, 'Screensize'));
% data_dfof = cat(3,data_dfof_targ_noadapt, data_dfof_targ_fake, data_dfof_dir_all, data_dfof_targ); % concat adapter resp, baseline2, targ resp
data_dfof = cat(3, data_dfof_targ_fake, data_dfof_dir_all); % did not concat data_dfof_targ_noadapt and data_dfof_targ due to higher noise
% current data_dfof contains: 
    % R1 (adapter response - baseline 1), 
    % decayed R1 (baseline 2 - baseline 1), 
    % and fake R2 (target response - baseline 1), which includes ad-tg and noad-tg

myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof, myfilter), [], 3);
% data_dfof_perc = prctile(imfilter(data_dfof, myfilter), 90, 3); % take top percentile due to noisy max tif?

figure;
imagesc(data_dfof_max)
title('data dfof max')
% figure;
% imagesc(data_dfof_perc)
% title('data dfof perc')

data_dfof = cat(3, data_dfof, data_dfof_max);
% data_dfof_max = data_dfof_perc;
    
otherwise
    disp('no such stim. stim type should be bunny, mix, grat6 or grating')
end

%% export tif for cellpose

disp('saving dfof-max-gauss tif & mat')
save data_dfof.mat data_dfof data_dfof_max

save_mat_as_tif(data_dfof_max);

disp(['mouse ', num2str(mouse), ' date ', num2str(date), ' ', run_str ...
    ' - finished data_reg & cellpose tif'])

end