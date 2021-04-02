%% get path names
clear all
clc

date = '200803';
mouse = 'i1322';
ImgFolder = strvcat('002');
time = strvcat('1140');
frame_rate = 30;

doFromRef = 0;
ref = strvcat('001'); % what is ref?
nrun = size(ImgFolder,1);
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
    nframes = max([temp(irun).counterValues{end}(end) info.config.frames]);

    
    fprintf(['Reading run ' num2str(irun) ', consisting of ' num2str(min(nframes)) ' frames \r\n'])
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

%% Choose register interval
nep = floor(size(data,3)./10000);
[n, n2] = subplotn(nep);
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

select = 4;
start_idx = select * 10000 + 1;
stop_idx = select * 10000 + 500;
data_avg = mean(data(:,:,start_idx:stop_idx),3);

%% Register data

% if exist(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
%     load(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
%     save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
%     [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
%     clear out outs
% elseif doFromRef
%     ref_str = ['runs-' ref];
%     if size(ref,1)>1
%         ref_str = [ref_str '-' ref(size(ref,1),:)];
%     end
%     load(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
%     [out, data_reg] = stackRegister(data,data_avg);
%     mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
%     save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
%     %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
%     %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
%     save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
% else
    tic; [out, data_reg] = stackRegister(data,data_avg); toc;
%     mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
%     save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
%     save(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
% end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
set(gcf, 'Position', get(0, 'Screensize'));
% print(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
set(gcf, 'Position', get(0, 'Screensize'));
% print(fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% find activated cells

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

%% determine ca signal latency (around 8 frames for this session)

tc_screen = mean(mean(data_reg,1),2);
tc_screen = squeeze(tc_screen);
all_trial_len = sz(3) - cStart(1);

data_trial = zeros(200, nTrials); % take 1-200 frame of every trial
data_trial_real = zeros(max(trial_len), nTrials);
whos tc_screen

trial_len = diff(cStart);
for it = 1:(nTrials-1)
    start_id = cStimOn(it);
    data_trial(:,it) = tc_screen(start_id : start_id + 200 - 1);
    data_trial_real(:,it) = [tc_screen(start_id : start_id + trial_len(it) - 1); NaN(max(trial_len) - trial_len(it), 1)];
end

plot(mean(data_trial, 2))
set(gcf, 'Position', get(0, 'Screensize'));
analysis_dir = fullfile(LL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);

data_trial_zoom_in = nanmean(data_trial_real, 2); plot(data_trial_zoom_in(1:50)); grid on; grid minor
set(gcf, 'Position', get(0, 'Screensize'));

%% fake movie
% trial type = [noad, 750, 250] x [ori 1-8]
% take original OR dfof movie to average over reps of same trial type

save_flag = 0; % toggle this to save/skip all .mat creation below
clear id_ad id_noad id_isi2 id_isi3 id_ori
clear frame_rate range_base range_resp ncell ntrial trial_len_max nisi nori ori_list
global id_ad id_noad id_isi2 id_isi3 id_ori % declare all global var for single dataset
global frame_rate range_base range_resp ncell ntrial trial_len_max nisi nori ori_list

imouse = mouse; area = 'V1';
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');
[~, frame_rate, input_behav, info, ~] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area);

% % data_reg = ny * nx * nframe
% npSub_tc = reshape(data_reg, [size(data_reg,1)*size(data_reg,2), size(data_reg,3)]);
% clear data_reg
% npSub_tc = npSub_tc'; % fake npSub_tc = nframe x npixel 
% save('npSub_tc_fake.mat', 'npSub_tc', '-v7.3') % force save >2GB .mat
tic; load npSub_tc_fake.mat; toc

%% params & indexing trials
% index by adapter contrast, target ori, isi

ntrial = input_behav.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps 
% final trial discarded bc too few frames
[nframe, ncell] = size(npSub_tc);

contrast_ad = celleqel2mat_padded(input_behav.tBaseGratingContrast); 
id_noad = find(contrast_ad == 0); id_ad = find(contrast_ad == 1); 
id_noad(id_noad > ntrial) = []; id_ad(id_ad > ntrial) = []; 

frame_ad = double(cell2mat(input_behav.cStimOn)); frame_ad_off = double(cell2mat(input_behav.cStimOff));
frame_tg = celleqel2mat_padded(input_behav.cTargetOn); frame_tg = double(frame_tg);
isi_seq = frame_tg - frame_ad_off; 
nisi = length(unique(frame_tg - frame_ad));
id_750 = find(isi_seq > mean(isi_seq)); id_250 = find(isi_seq < mean(isi_seq)); 
id_750(id_750 > ntrial) = []; id_250(id_250 > ntrial) = []; 
id_ad750 = intersect(id_ad, id_750); id_ad250 = intersect(id_ad, id_250);
id_isi2 = {id_ad750, id_ad250}; 
id_isi3 = {id_noad, id_ad750, id_ad250};

trial_len_max = max(unique(diff(frame_ad)));

ori_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg); ori_seq(ori_seq == 180) = 0;
ori_seq(end) = [];
ori_list = unique(ori_seq); 
nori = length(ori_list); id_ori = cell(nori, 1);
for iori  = 1 : nori
    id_ori{iori} = find(ori_seq == ori_list(iori)); 
end

%% 
tic
tc_align_ad = align_tc(frame_ad, npSub_tc);
toc

%%
tic
save('tc_align_ad.mat', 'tc_align_ad', '-v7.3') % force save >2GB .mat
toc

range_base = 1:3; range_resp = 9:12;

tic
[trace_avg, trace_sem] = trace_grand_avg(tc_align_ad, save_flag);
toc