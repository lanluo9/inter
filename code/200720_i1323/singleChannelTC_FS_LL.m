%% get path names
date = '200720';
ImgFolder = strvcat('003');
time = strvcat('1133');
mouse = 'i1323';
doFromRef = 0;
ref = strvcat('002'); % what is ref?
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LL_base '\Data\2P_images\' mouse '\' date '\' ImgFolder(irun,:)];
%     CD = [LL_base '\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];

    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
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
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc

% %% For behavior experiments
% % Plot outcome by trial number
% SIx = strcmp(input.trialOutcomeCell, 'success');
% MIx = strcmp(input.trialOutcomeCell, 'ignore');
% 
% figure;
% plot(smooth(SIx,10));
% hold on
% plot(smooth(MIx,10));
% 
% % Crop data and input struct
% input = trialChopper(input,[1 200]);
% data = data(:,:,input.counterValues{1}(1):input.counterValues{end}(end));

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

data_avg = mean(data(:,:,70001:70500),3);

%% Register data

if exist(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    tic; [out, data_reg] = stackRegister(data,data_avg); toc;
    mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'data_reg')
    save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
print(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% find activated cells

% tCyc = cell2mat(input.tCyclesOn);
cStart = cell2mat(input.cStimOn); % same as cStimOn
cStimOn = cell2mat(input.cStimOn);
cStimOff = cell2mat(input.cStimOff);
cTarget = celleqel2mat_padded(input.cTargetOn); cTarget = int64(cTarget);
nTrials = input.trialsSinceReset; % 464 = 32 types * 14.5 reps
sz = size(data_reg); % [y pixel * x pixel * nframe]

data_f = zeros(sz(1),sz(2),nTrials);
data_adapter = zeros(sz(1),sz(2),nTrials);
data_f2 = zeros(sz(1),sz(2),nTrials);
data_targ = zeros(sz(1),sz(2),nTrials);

%% trial design
unique(cStimOff - cStimOn) % adapter = 3-4 frame. rounded from 100 ms?
unique(cTarget - cStimOff) % ISI = 8 or 22-23 frame. rounded from 250? or 750 ms
unique(cTarget - cStimOn)
% target: 3-4 frame too

% index2 = structfun(@(x) any(contains(x, 'Targ')), input);  sum(index2)
tt = fieldnames(input)
index = cellfun(@(x) any(contains(x, 'targ')),tt); sum(index)
id = find(index > 0);
for i = 1 : length(id)
    fprintf(['input.', tt{id(i)}])
    disp(' ')
end

input.cTargetOn 
% input.tSoundTargetAmplitude 
% input.doLinearTargetSpacing 
% input.soundTargetAmplitude 
% input.soundTargetStepsPerOctave 
input.block2TargetOnTimeMs % 100 ms -> thus also 3-4 frames
% input.block2SoundTargetAmplitude 
% input.block2SoundTargetStepsPerOctave 
input.doCustomTargetTime % 0
input.targetStimOnMs
input.targetOnTimeMs % 100 ms

% %%
% 
trial_len = diff(cStart);
unique(trial_len)
histogram(trial_len)
% 
% pretend_len = 207;
% tt = data_reg(:,:, 1 : sz(3) - mod(sz(3), pretend_len));
% t = sz(3) - mod(sz(3), pretend_len);
% tt = reshape(tt, [sz(1)*sz(2), pretend_len, t/pretend_len]);
% tt_avg = mean(tt, 3);
% 
% imagesc(tt_avg(1:1000, :))

%% determine ca signal latency (around 8 frames for this session)

tc_screen = mean(mean(data_reg,1),2);
tc_screen = squeeze(tc_screen);
all_trial_len = sz(3) - cStart(1);

data_trial = zeros(200, nTrials); % take 1-200 frame of every trial
data_trial_real = zeros(max(trial_len), nTrials);
whos tc_screen

for it = 1:(nTrials-1)
    start_id = cStimOn(it);
    data_trial(:,it) = tc_screen(start_id : start_id + 200 - 1);
    data_trial_real(:,it) = [tc_screen(start_id : start_id + trial_len(it) - 1); NaN(max(trial_len) - trial_len(it), 1)];
end

plot(mean(data_trial, 2))
data_trial_zoom_in = nanmean(data_trial_real, 2); plot(data_trial_zoom_in(1:50)); grid on; grid minor

% %%
% 
% 
% tt = tc_screen(cStart(1) : all_trial_len - mod(all_trial_len-cStart(1)+1, pretend_len));
% temp = reshape(tt, [pretend_len, length(tt)/pretend_len]);
% figure
% plot(mean(temp, 2))

%% calculate response
% see data_trial_zoom_in. ca signal latency = 8 frames & adapter|target = 3|4 frames
% count from frame #1
% data_adapter = frame #8-11
% data_f2 (baseline after adaptation) = frame #14-16

assert(length(cTarget) == nTrials && length(cStart) == nTrials && cTarget(nTrials)+3 < sz(3))
for itrial = 1:nTrials
%     if ~isnan(cStart(itrial))
        data_f(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)-10:cStart(itrial)-1),3);
        data_adapter(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+7:cStimOn(itrial)+10),3);
        
%         if cStimOn(itrial) >= cStart(itrial) 
        data_f2(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+13:cStimOn(itrial)+15),3);
%         else
%             data_base2(:,:,itrial) = nan(sz(1),sz(2));
%         end
%     else
%         data_bef(:,:,itrial) = nan(sz(1),sz(2));
%         data_adapter(:,:,itrial) = nan(sz(1),sz(2));
% %         data_base2(:,:,itrial) = nan(sz(1),sz(2));
%     end
    
%     if ~isnan(cTarget(itrial))
%         if cTarget(itrial)+3 < sz(3)
        data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+7:cTarget(itrial)+10),3);
%         else
%             data_targ(:,:,itrial) = nan(sz(1),sz(2));
%         end
%     else
%         data_targ(:,:,itrial) = nan(sz(1),sz(2));
%     end
end

data_adapter_dfof = (data_adapter-data_f)./data_f;
data_base2_dfof = (data_f2-data_f)./data_f;
data_targ_dfof = (data_targ-data_f2)./data_f2; 
data_targ_dfof_fake = (data_targ-data_f)./data_f; 

% %% find direction
% 
% tt = fieldnames(input)
% index = cellfun(@(x) any(contains(x, 'Dir')),tt); sum(index)
% id = find(index > 0);
% for i = 1 : length(id)
%     fprintf(['input.', tt{id(i)}])
%     disp(' ')
% end
% 
% % input.tGratingMaxDirectionStepDeg 
% input.tGratingDirectionStepsPerOctave % ??
% % input.tBaseGratingDirectionDeg 
% input.tGratingDirectionDeg % == input.gratingDirectionDeg
% % input.tCatchGratingDirectionDeg 
% % input.baseGratingDirectionDeg 
% % input.baseGratingDirectionStepDeg 
% % input.baseGratingDirectionStepN 
% % input.gratingMaxDirectionStepDeg 
% % input.gratingDirectionStepsPerOctave 
% % input.block2BaseGratingDirectionDeg % 0
% % input.block2GratingMaxDirectionStepDeg 
% % input.block2GratingDirectionStepsPerOctave 

%% plot response
targCon = celleqel2mat_padded(input.tGratingContrast);
unique(targCon) % target contrast 1
if input.doRandCon
    adapterCon = ones(size(targCon));
else
    adapterCon = celleqel2mat_padded(input.tBaseGratingContrast);
end
unique(adapterCon) % adapter contrast 0 or 1
ind_con = intersect(find(targCon == 1),find(adapterCon == 0));

adapterDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(adapterDir);
ndir = length(dirs);
targetDelta = celleqel2mat_padded(input.tGratingDirectionDeg);
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
imagesc(data_dfof2_dir) % but baseline2 shows bright cells, almost same as adapter resp
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
% data_dfof = cat(3,data_dfof_dir_all, data_dfof_targ); % concat adapter resp, baseline2, targ resp

data_dfof_targ_noadapt = zeros(sz(1),sz(2),nDelta);
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
data_dfof = cat(3,data_dfof_targ_fake, data_dfof_targ_noadapt, data_dfof_dir_all, data_dfof_targ); % concat adapter resp, baseline2, targ resp
% noadapt should be upfront!

myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof, myfilter),[],3);
figure;
imagesc(data_dfof_max)
title('data dfof max')

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;
% mask_data = data_dfof_targ_fake;

for iStim = 1:size(data_dfof,3)
% for iStim = 1:size(data_dfof_targ_fake,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    
    fprintf('%d out of %d',iStim, size(data_dfof,3));
    bwout = imCellEditInteractiveLG_LL(mask_data_temp);
%     bwout = imCellEditInteractive(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end

mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)
set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
saveas(gcf, ['mask_cell_addfake.jpg'])

% bwout = imCellEditInteractive(data_dfof_max);
% mask_cell = bwlabel(bwout);

%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell_addfake.mat']), 'data_dfof', 'mask_cell', 'mask_np')
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5;
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

save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs_addfake.mat']), 'data_tc', 'np_tc', 'npSub_tc')
clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

% %% FS cycle analysis
% 
% if iscell(input.nFramesOn)
%     nOn = input.nFramesOn{1};
% else
%     nOn = input.nFramesOn;
% end
% prewin_frames = 30;
% postwin_frames = 90;
% tCyc = cell2mat(input.tCyclesOn);
% cStart = celleqel2mat_padded(input.cFirstStim);
% cTarget = celleqel2mat_padded(input.cTargetOn);
% nTrials = length(tCyc);
% nCells = size(npSub_tc,2);
% maxCyc = max(tCyc,[],2);
% data_trial = nan(prewin_frames+postwin_frames,nCells,maxCyc+1,nTrials);
% 
% tFramesOff = nan(nTrials,maxCyc);
% SIx = strcmp(input.trialOutcomeCell, 'success');
% MIx = strcmp(input.trialOutcomeCell, 'ignore');
% FIx = strcmp(input.trialOutcomeCell, 'failure');
% nCyc = tCyc;
% nCyc([find(MIx) find(SIx)]) = tCyc([find(MIx) find(SIx)])+1;
% for itrial = 1:nTrials
%     if isfield(input, 'tFramesOff')
%         if length(input.tFramesOff{itrial}>0)
%             tempFramesOff = input.tFramesOff{itrial};
%         else
%             tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
%             input.tFramesOff{itrial} = tempFramesOff;
%         end
%     else
%         if iscell(input.nFramesOff)
%             tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
%         else
%             tempFramesOff = input.nFramesOff.*(ones(1,tCyc(itrial)));
%         end
%     end
% 
%     tFramesOff(itrial,1:tCyc(itrial)) = tempFramesOff(1:tCyc(itrial));
%     if ~isnan(cStart(itrial))
%         for icyc = 1:nCyc(itrial)
%             if icyc > 1
%                 cyc_add = ((icyc-1)*nOn)+sum(tempFramesOff(1:icyc-1));
%             else
%                 cyc_add = 0;
%             end
%             if cStart(itrial)+postwin_frames-1+cyc_add <= size(npSub_tc,1)
%                 data_trial(:,:,icyc,itrial) = npSub_tc(cStart(itrial)-prewin_frames+cyc_add:cStart(itrial)+postwin_frames+cyc_add-1,:);
%             else
%                 data_trial(:,:,icyc,itrial) = NaN(prewin_frames+postwin_frames,nCells);
%             end 
%         end
%     else
%         data_trial(:,:,icyc,itrial) = NaN(prewin_frames+postwin_frames,nCells);
%     end
% end
% data_f = nanmean(data_trial(1:prewin_frames,:,1,:),1);
% data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);
% 
% targCon = celleqel2mat_padded(input.tGratingContrast);
% if isfield(input,'doRandCon') & input.doRandCon
% 	adapterCon = nan(maxCyc,nTrials);
%     for itrial = 1:nTrials
%         adapterCon(:,itrial) = input.tBaseGratingContrast{itrial}(1:tCyc(itrial));
%     end
%     ind_con = [];
% else
%     adapterCon = celleqel2mat_padded(input.tBaseGratingContrast);
%     ind_con = intersect(find(targCon == 1),find(adapterCon == 0));
% end
% adapterDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
% dirs = unique(adapterDir);
% ndir = length(dirs);
% tGratingDir = round(double(celleqel2mat_padded(input.tGratingDirectionDeg)),0);
% if sum(tGratingDir-adapterDir) == 0
%     targetDelta = tGratingDir-adapterDir;
% else
%     targetDelta = tGratingDir;
% end
% deltas = unique(targetDelta);
% nDelta = length(deltas);
% offs = unique(tFramesOff(:,1));
% noff = length(offs);
% frameRateHz = input.frameRateHz;
% 
% base_win =33:35;
% resp_win =39:41; 
% 
% % figure;
% % if nCells<25
% %     ii = nCells;
% % else
% %     ii = 25;
% % end
% % for i = 1:ii
% %     subplot(5,5,i)
% % if length(ind_con)>10
% %     plot(squeeze(nanmean(mean(data_dfof(20:50,i,2,ind_con),2),4)))
% % elseif noff>1
% %     ind = find(tFramesOff(:,1) == offs(noff));
% %     plot(squeeze(nanmean(mean(data_dfof(20:50,i,1,:),2),4)))
% % else
% %     plot(squeeze(nanmean(mean(data_dfof(20:50,i,1,:),2),4)))
% % end
% % vline(base_win-19)
% % vline(resp_win-19)
% % end
% 
% figure;
% subplot(2,1,1)
% plot(squeeze(nanmean(mean(data_dfof(:,:,1,:),2),4)));
% vline(base_win,'k:')
% vline(resp_win,'r:')
% title('Baseline')
% subplot(2,1,2)
% sz = size(data_dfof);
% data_targ = zeros(sz(1),sz(2),length([find(SIx)]));
% for itrial = 1:sz(4);
%     %if find([find(SIx)] == itrial)
%         data_targ(:,:,itrial) = data_dfof(:,:,nCyc(itrial),itrial);
%     %end
% end
% plot(squeeze(nanmean(mean(data_targ,2),3)));
% title('Target')
% vline(base_win,'k:')
% vline(resp_win,'r:')
% 
% %%
% 
% save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof', 'prewin_frames', 'postwin_frames', 'base_win','resp_win')
% save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'prewin_frames','baseDir', 'dirs', 'ndir', 'tFramesOff', 'offs', 'noff', 'baseCon', 'ind_con', 'tGratingDir', 'targetDelta', 'deltas', 'nDelta', 'tCyc', 'nCyc', 'maxCyc', 'nCells', 'frameRateHz', 'nTrials', 'SIx', 'MIx', 'FIx', 'cTarget', 'cStart', 'base_win','resp_win')
