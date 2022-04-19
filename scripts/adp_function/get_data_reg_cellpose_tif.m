function [data_reg, LL_base, date, imouse, run_str] = get_data_reg_cellpose_tif(arg_mouse, arg_date, arg_ImgFolder)

disp('start running get_data_reg_cellpose_tif.m')
mouse = arg_mouse
imouse = ['i', num2str(mouse)];
date = num2str(arg_date)
ImgFolder = arg_ImgFolder

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images';
xls_dir = fullfile(fn_base, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx'); clear dataset_meta
dataset_meta = readtable(xls_file.name); 
idx = find(all(ismember(dataset_meta.(1),[ImgFolder,'_000_000']),2));
time = num2str(dataset_meta.(8)(idx));
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);

% load and register
disp('start loading sbx')
tic
data = [];
clear temp
trial_n = [];
offset = 0;

for irun = 1:nrun
    LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
    CD = [LL_base '\Data\2P_images\' imouse '\' date '\' ImgFolder(irun,:)];
    cd(CD);
    
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' imouse '-' date '-' time(irun,:) '.mat'];
    behav_mat = load(fName);
    temp(irun) = behav_mat.input;
    behav_input = behav_mat.input;
    nframes = max([temp(irun).counterValues{end}(end) info.config.frames]);
    fprintf(['Reading run ' num2str(irun) ', consisting of ' num2str(min(nframes)) ' frames \r\n'])

    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    offset = offset+min(nframes);
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
% behav_input = concatenateDataBlocks(temp);
clear data_temp temp
toc

% Choose register interval
nep = floor(size(data,3)./10000);
[n, n2] = subplotn(nep);
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep
    subplot(n,n2,i); 
    imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); 
    title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); 
end
disp('start registration. using middle section as motion correct ref')
select = 4
start_idx = select * 10000 + 1;
stop_idx = select * 10000 + 500;
data_avg = mean(data(:,:,start_idx:stop_idx),3);

% Register data
if exist(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'behav_input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
else
    tic; [out, data_reg] = stackRegister(data,data_avg); toc;
%     mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'behav_input')
end
clear data out
disp('registration finished')

% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

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

assert(behav_input.doRandStimOne == 1 & behav_input.doSameStims == 1)
adapter_id = cell2mat(behav_input.tstimOne);
adapter_list = unique(adapter_id);
n_adapter = length(adapter_list);
target_id = cell2mat(behav_input.tstimTwo);
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
close all

%% export tif for cellpose

stim_resp_gauss = data_dfof_max'; % gauss smooth each stim resp, take max
disp('saving dfof-max-gauss tif & mat')
save data_dfof.mat data_dfof data_dfof_max

tif_file = ['cellpose_stim_resp_gauss.tif'];
fTIF = Fast_BigTiff_Write(tif_file,1,0);
tic
msg = 0;
N = 1;
B = numel(stim_resp_gauss)*2;
for ct = 1:N
    fprintf(1,repmat('\b',[1,msg]));
    msg = fprintf(1,'%.0f/%.0f',ct,N);
    fTIF.WriteIMG(stim_resp_gauss(:,:,ct));
end
fTIF.close;
fprintf(1,repmat('\b',[1,msg]));
t=toc
fprintf(1,'\nWrite %.0f bytes in %.0f mins \n',B*N,t/60);
fprintf(1,'Write speed: %.0f MB/s \n',(B*N)/(2^20*t));

disp(['mouse', num2str(mouse), 'date', num2str(date), ...
    'finished data_reg & cellpose tif'])

end