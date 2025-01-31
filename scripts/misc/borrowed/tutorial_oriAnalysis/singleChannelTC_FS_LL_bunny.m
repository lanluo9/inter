%% get path names
clear all
clc

%%
database_path = 'Z:\All_Staff\home\lan\Data\2P_images\';
master_xls = [database_path, 'mat_inter\adp_dataset_master.xlsx'];
dataset_meta = readtable(master_xls);
dataset_now = dataset_meta(ismember(dataset_meta.paradigm, ...
    'bunnytop high res'),:); % high lum-contrast
dataset_now = dataset_now(dataset_now.mouse == 1369, :);
dataset_now = dataset_now(ismember(dataset_now.area, ...
    'V1'),:);
dataset_now

iset = 1
mouse = dataset_now.mouse(iset)
imouse = ['i', num2str(mouse)];
date = num2str(dataset_now.date(iset))
ImgFolder = dataset_now.num(iset); ImgFolder = ImgFolder{1}

%%
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images';
xls_dir = fullfile(fn_base, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx'); clear dataset_meta
dataset_meta = readtable(xls_file.name); 
idx = find(all(ismember(dataset_meta.(1),[ImgFolder,'_000_000']),2));
time = num2str(dataset_meta.(8)(idx));
frame_rate = 30;

doFromRef = 0;
ref = strvcat('001'); % what is ref?
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
    CD = [LL_base '\Data\2P_images\' imouse '\' date '\' ImgFolder(irun,:)];
%     CD = [LL_base '\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' imouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = max([temp(irun).counterValues{end}(end) info.config.frames]);

    
    fprintf(['Reading run ' num2str(irun) ', consisting of ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
%     data_temp = sbxread(imgMatFile(1,1:11),0,100000);
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cLeverUp') 
        if irun>1
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                %temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                temp(irun).cLeverDown{itrial} = temp(irun).cLeverDown{itrial}+offset;
                temp(irun).cFirstStim{itrial} = temp(irun).cFirstStim{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                if ~isempty(temp(irun).cLeverUp{itrial})
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial}+offset;
                else
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial};
                end
                if ~isempty(temp(irun).cTargetOn{itrial})
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial}+offset;
                else
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial};
                end
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
% input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc

% %% Choose register interval
nep = floor(size(data,3)./10000);
[n, n2] = subplotn(nep);
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

%%
select = 3
start_idx = select * 10000 + 1;
stop_idx = select * 10000 + 500;
data_avg = mean(data(:,:,start_idx:stop_idx),3);

% %% Register data

if exist(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]))
    load(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' ref_str], [date '_' imouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'input')
else
    tic; [out, data_reg] = stackRegister(data,data_avg); toc;
    mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'input')
end
clear data out

% %% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% find activated cells

% tCyc = cell2mat(input.tCyclesOn);
cStart = cell2mat(input.cStimOneOn); % same as cStimOn
cStimOn = cell2mat(input.cStimOneOn);
cStimOff = cell2mat(input.cStimOneOff);
% cTarget = celleqel2mat_padded(input.cStimTwoOn); cTarget = int64(cTarget);
cTarget = cell2mat(input.cStimTwoOn);
try nTrials = input.trialsSinceReset % 464 = 32 types * 14.5 reps
catch
    nTrials = input.trialSinceReset
end
sz = size(data_reg); % [y pixel * x pixel * nframe]

data_f = zeros(sz(1),sz(2),nTrials);
data_adapter = zeros(sz(1),sz(2),nTrials);
data_f2 = zeros(sz(1),sz(2),nTrials);
data_targ = zeros(sz(1),sz(2),nTrials);

% %% trial design
unique(cStimOff - cStimOn) % adapter = 3-4 frame. rounded from 100 ms?
unique(cTarget - cStimOff) % ISI = 8 or 22-23 frame. rounded from 250? or 750 ms
unique(cTarget - cStimOn) % target: 3-4 frame too

trial_len = diff(cStart); unique(trial_len)
% histogram(trial_len)

% %% determine ca signal latency (around 8 frames for this session)

tc_screen = mean(mean(data_reg,1),2);
tc_screen = squeeze(tc_screen);
all_trial_len = sz(3) - cStart(1);

data_trial = zeros(min(unique(trial_len)), nTrials); % take 1-200 frame of every trial
data_trial_real = zeros(max(trial_len), nTrials);
whos tc_screen

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
disp('now manually update ca latency and adapted baseline window in next section')
disp(' ')

%% calculate response
% see data_trial_zoom_in
% if ca signal latency = +7 frames 
% adapter|target = 3|4 frames
% count from frame #1
% data_adapter = frame #8-11
% data_f2 (baseline after adaptation) = frame #14-16

% ca_latency = 5 or 8;
ca_latency = 13; % = x-1. stim onset frame 1 -> signal received frame x
window_len = 3; % 2-3

% ca_latency = 4;
% window_len = 1;

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

assert(input.doRandStimOne == 1 & input.doSameStims == 1)
adapter_id = cell2mat(input.tstimOne);
adapter_list = unique(adapter_id);
n_adapter = length(adapter_list);
target_id = cell2mat(input.tstimTwo);
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
data_dfof_max = max(imfilter(data_dfof, myfilter),[],3);

% figure
% imagesc(max(data_dfof_ad, [], 3)) 
% title('data dfof ad max')
% set(gcf, 'Position', get(0, 'Screensize'));
% figure
% imagesc(max(data_dfof_targ, [], 3)) 
% title('data dfof tg max')
% set(gcf, 'Position', get(0, 'Screensize'));
% figure; imagesc(data_dfof_max)
% title('data dfof max')

data_dfof = cat(3,data_dfof, data_dfof_max); % adapter, targ, targ_fake, gaussian filter max proj

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    
    fprintf('\n %d out of %d \n',iStim, size(data_dfof,3));
    bwout = imCellEditInteractiveLG_LL(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end

mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['mask_cell_addfake.jpg'])
disp('mask cell img saved')

%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_mask_cell_addfake.mat']), 'data_dfof', 'mask_cell', 'mask_np')
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5; % down sampling
sz = size(data_reg);
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2)

np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s /n']) 
     disp(' ')
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew, ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_TCs_addfake.mat']), 'data_tc', 'np_tc', 'npSub_tc')
clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
disp('TC extraction complete')
