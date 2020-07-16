%% 

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey'); 
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); % mwork = behavior data
tc_fn = fullfile(lg_fn, 'Analysis\2P');
fnout = fullfile(lg_fn, 'Analysis\2P\test'); % output in analysis folder

date = '200118';
ImgFolder = '002';
time = strvcat('1508'); % catenate multiple time vertically to get a matrix of file/folder name
mouse = 'i1312';
frame_rate = 15.5;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%% load data

fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName); % load behavior data, aka "input"

CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata, aka "info" % check content

tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs.mat']);

%% use npSub_tc to conduct further analysis

nOn = input.nScansOn; % behavior data "input"
nOff = input.nScansOff;
trial_len = nOn + nOff;
ntrials = size(input.tGratingDirectionDeg,2);
size(npSub_tc) % nframe * ncell
ncell = size(npSub_tc, 2);

tc_trials = reshape(npSub_tc,[nOn+nOff, ntrials, size(npSub_tc,2)]); 
tc_trials = permute(tc_trials, [3,2,1]);
size(tc_trials) % ncell * ntrial * trial_len

Dir = celleqel2mat_padded(input.tGratingDirectionDeg); 
convert_idx = Dir>=180;
Ori = Dir;
Ori(convert_idx) = Ori(convert_idx) - 180; 
Ori_list = unique(Ori);
nOri = length(Ori_list);

% tt = mean(mean(tc_trials,1),2);
% plot(1:trial_len, squeeze(tt)) 
% trial = 0-60 off + 61-90 on. signal decays >2/3 after 0-30 off
% as short as 3 frames would suffice. now take 10 frames as window len:
win_len = 10;

%% cells sensitive to orientations
sig_ttest = pi * ones(ncell, nOri);
resp_avg = pi * ones(ncell, nOri);
resp_ste = pi * ones(ncell, nOri); % standard error 
dfof_avg = pi * ones(ncell, nOri); % dF/F
dfof_ste = pi * ones(ncell, nOri);

for iOri = 1 : nOri
    idx = find(Ori == Ori_list(iOri)); 
    ntrials_ori = length(idx);
    for icell = 1 : ncell
        base_win = squeeze(tc_trials(icell, idx, (nOff - win_len + 1):nOff));
        base_win = mean(base_win, 2); % avg over window -> [ntrial_ori, 1]
        resp_win = squeeze(tc_trials(icell, idx, (trial_len - win_len + 1):trial_len));
        resp_win = mean(resp_win, 2);
        
        sig_ttest(icell, iOri) = ttest(base_win, resp_win, 'alpha',0.05./(ntrials_ori));
        resp_avg(icell, iOri) = mean(resp_win);
        dfof_avg(icell, iOri) = mean( (resp_win - base_win) ./ mean(base_win) );
        resp_ste(icell, iOri) = std(resp_win) / sqrt(length(resp_win));
        dfof_ste(icell, iOri) = std( (resp_win - base_win) ./ mean(base_win) ) / sqrt(ntrials_ori);
    end
end

sum(sum(sig_ttest,2)>0) % ncells responsive to >= 1 orientation: 80/148

size(tc_trials) % ncell * ntrial * trial_len
base = mean(tc_trials(:,:, (nOff - win_len + 1):nOff), 3);
resp = mean(tc_trials(:,:, (trial_len - win_len + 1):trial_len), 3);
df = resp - base;

%% orientation tuning of indiv cell

for icell = 1 : 25 % ncell
    figure('units','normalized','outerposition',[0 0 1 1]);
    errorbar(Ori_list, dfof_avg(icell,:), dfof_ste(icell,:), 'LineWidth',1)
    hold on
    if sum(sig_ttest(icell,:)) > 0
        sig_idx = sig_ttest(icell,:) > 0;
        sig_ori = Ori_list(sig_idx);
        sig_star_height = dfof_avg(icell,sig_idx) + dfof_ste(icell,sig_idx) + 0.01;
        scatter(sig_ori, sig_star_height, '*', 'LineWidth',1)
    else
        sig_ori = [];
    end
    xlim([0-5, 180])
    title(['cell ', num2str(icell), ': sensitive to ' num2str(length(sig_ori)), ' orientations'])
    saveas(gcf, ['ori_tuning_', num2str(icell)], 'jpg')
    close
end
% cell 2 & 17?
% sig & dfof correct?


%%
