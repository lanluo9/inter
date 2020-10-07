%% set up 

close all
clear
clc
cd C:\Users\lan\Documents\repos\inter\mat

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806,...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};

%% iterate thru datasets

for iset = 1 : length(dataset_list.date)
iset
date = num2str(dataset_list.date(iset))
mouse = num2str(dataset_list.mouse(iset))
area = dataset_list.area{1,iset}

imouse = ['i', mouse];
xls_dir = fullfile(data_fn, imouse, date);
cd(xls_dir)
xls_file = dir('*.xlsx'); 
clear dataset_meta
dataset_meta = readtable(xls_file.name); 
time = dataset_meta.(8){end}
ImgFolder = dataset_meta.(1){end}(1:3);
frame_rate = str2num(dataset_meta.(5){end});

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' imouse];
datemouserun = [date '_' imouse '_' run_str];
areamousedate = [dataset_list.area{1,iset} '_' imouse '_' date];

fName = fullfile(mworks_fn, ['data-' imouse '-' date '-' time '.mat']);
load(fName); % load behavior data "input"
input_behav = input;
clear input

CD = fullfile(data_fn, imouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata "info"

tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs_addfake.mat']); % load time course including fake targ resp
% fix bug later: retinotopy folder in Analysis naming convention should adhere to tc folder

result_prefix = 'C:\Users\lan\Documents\repos\inter\mat\';
result_folder = fullfile(result_prefix, areamousedate);
if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder);

%% dataset params

ntrial = input_behav.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps % final trial discarded bc too few frames
[nframe, ncell] = size(npSub_tc); % nframe * ncell

cStimOn = double(cell2mat(input_behav.cStimOn)); cStimOff = double(cell2mat(input_behav.cStimOff));
cTarget = celleqel2mat_padded(input_behav.cTargetOn); cTarget = double(cTarget);

adapter_stim_len = min(unique(cStimOff - cStimOn)); % adapter = 3-4 frames, rounded from 100 ms
isi_len = unique(cTarget - cStimOff); % ISI = 8 or 22-23 frames, rounded up from 250 or 750 ms
targ_stim_len = floor(input_behav.targetOnTimeMs / frame_rate); % target = 3-4 frames, rounded from 100 ms
iti_len = unique(cStimOn(2:end) - (cTarget(1:end-1) + 3)); % ITI = 192-193-194 frames, 6.4-6.5 s?
trial_len_list = unique(diff(cStimOn)); % 6.9 or 7.4s

targCon = celleqel2mat_padded(input_behav.tGratingContrast); unique(targCon); % target contrast 1
adapterCon = celleqel2mat_padded(input_behav.tBaseGratingContrast); unique(adapterCon); % adapter contrast 0 or 1

adapterDir = celleqel2mat_padded(input_behav.tBaseGratingDirectionDeg);
dirs = unique(adapterDir); % adapter dir === 0
ndir = length(dirs);
delta_seq = celleqel2mat_padded(input_behav.tGratingDirectionDeg);
delta_list = unique(delta_seq); % target 8 dir (actually ori): 22.5-180. equivalent to diff from adapter
% fix bug later: need to convert 180 to 0
ndelta = length(delta_list); 

%% index of trials

isi_list = cTarget - cStimOff; 
ngap = length(unique(cTarget - cStimOn));

id_750 = find(isi_list > mean(isi_list)); id_750(id_750 > ntrial) = [];
id_250 = find(isi_list < mean(isi_list)); id_250(id_250 > ntrial) = [];
id_gaps = {id_750, id_250}; 
id_noad = intersect(find(targCon == 1),find(adapterCon == 0)); id_noad(id_noad > ntrial) = [];
id_ad = intersect(find(targCon == 1),find(adapterCon == 1)); id_ad(id_ad > ntrial) = [];

%% align tc by adapter onset or targ onset. normalize aligned tc by 1-sec-long "trial baseline"

trial_len_ad = diff(cStimOn);
tc_trial_align_ad = zeros(ncell, ntrial, max(trial_len_ad));
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cStimOn(itrial);
        tc_trial_align_ad(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + trial_len_ad(itrial) - 1); ...
            NaN(max(trial_len_ad) - trial_len_ad(itrial), 1)];
    end
end

trial_len_targ = diff(cTarget);
tc_trial_align_targ = zeros(ncell, ntrial, max(trial_len_targ));
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cTarget(itrial);
        tc_trial_align_targ(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + trial_len_targ(itrial) - 1); ...
            NaN(max(trial_len_targ) - trial_len_targ(itrial), 1)];
    end
end

trial_base_len = frame_rate * 1; % 30 frame/sec * 1 sec
tc_trial_base = zeros(ncell, ntrial, trial_base_len);
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cStimOn(itrial) - trial_base_len;
        tc_trial_base(icell, itrial, :) = [npSub_tc_cell(cStimOn(itrial) - trial_base_len : cStimOn(itrial) - 1)];
    end
end
t = squeeze(nanmean(squeeze(tc_trial_base(:,:,:)), 1)); 
t_base = squeeze(nanmean(t(:,:), 1)); 
plot(t_base, 'k');

result_sub = fullfile(result_folder, 'pre-processing');
if ~exist(result_sub); mkdir(result_sub); end
cd(result_sub)
% saveas(gcf, ['trial base'], 'jpg'); close 

tc_trial_base_avg = mean(tc_trial_base, 3);
for icell = 1:ncell
    for itrial = 1:ntrial
        tc_trial_align_ad(icell, itrial, :) = tc_trial_align_ad(icell, itrial, :) ./ tc_trial_base_avg(icell, itrial, :) -1;
        tc_trial_align_targ(icell, itrial, :) = tc_trial_align_targ(icell, itrial, :) ./ tc_trial_base_avg(icell, itrial, :) -1;
    end
end

%% find base & resp window in aligned tc
% range = 223; 
range = 50;
t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,:,:)), 1));
t_ad = squeeze(nanmean(t(:,:), 1)); 
t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,:,:)), 1)); 
t_tg = squeeze(nanmean(t(:,:), 1)); 

plot(t_ad(1:range), 'r'); hold on; plot(t_tg(1:range), 'b'); 
grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
legend('ad align', 'targ align')
saveas(gcf, ['aligned_tc_zoomin'], 'jpg')
close

% t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,:,:)), 1));
% t_ad = squeeze(nanmean(t(:,:), 1)); 
% t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_750,:)), 1)); 
% t_tg_750 = squeeze(nanmean(t(:,:), 1)); 
% t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_250,:)), 1)); 
% t_tg_250 = squeeze(nanmean(t(:,:), 1)); 
% 
% plot(t_ad(1:range), 'r'); hold on; plot(t_tg_750(1:range), 'b'); plot(t_tg_250(1:range), 'g'); 
% grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
% legend('ad align', 'targ align 750', 'targ align 250')

% % t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,id_750,:)), 1));
% % t_ad_750 = squeeze(nanmean(t(:,:), 1)); 
% % t = squeeze(nanmean(squeeze(tc_trial_align_ad(:,id_250,:)), 1));
% % t_ad_250 = squeeze(nanmean(t(:,:), 1)); 
% % 
% % t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_750,:)), 1)); 
% % t_tg_750 = squeeze(nanmean(t(:,:), 1)); 
% % t = squeeze(nanmean(squeeze(tc_trial_align_targ(:,id_250,:)), 1)); 
% % t_tg_250 = squeeze(nanmean(t(:,:), 1)); 
% % 
% % plot(t_ad_750(1:range), 'r'); hold on; plot(t_ad_250(1:range), 'r--'); 
% % plot(t_tg_750(1:range), 'b'); plot(t_tg_250(1:range), 'g'); 
% % grid on; grid minor; set(gcf, 'Position', get(0, 'Screensize'));
% % legend('ad align 750', 'ad align 250', 'targ align 750', 'targ align 250')


range_base = [1:3]; 
range_resp = [9:12];
% prompt = 'base window = 1:3. what is resp window? ';
% range_resp = input(prompt)
% close

%% find vis-driven & good-fit

% vis-driven by noad targ
sig_ttest_noad = pi * ones(ncell, ndelta); p_ttest_noad = pi * ones(ncell, ndelta);
base_avg_noad = pi * ones(ncell, ndelta);
resp_avg_noad = pi * ones(ncell, ndelta); resp_ste_noad = pi * ones(ncell, ndelta); % standard error 
dfof_avg_noad = pi * ones(ncell, ndelta); dfof_ste_noad = pi * ones(ncell, ndelta); % dF/F

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
for icell = 1 : ncell
    base_win = []; resp_win = [];
    
    for igap =  1 : ngap % order: 750, 250
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_noad); % use only no-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 
        
        base_win = [base_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_resp)),2)];
    end
        
    [sig_ttest_noad(icell, idelta), p_ttest_noad(icell, idelta)] = ttest(base_win, resp_win,...
            'alpha',0.01./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg_noad(icell, idelta) = mean(base_win); % avg over trials of same ori
    resp_avg_noad(icell, idelta) = mean(resp_win);
    resp_ste_noad(icell, idelta) = std(resp_win) / sqrt(length(resp_win));

    dfof_avg_noad(icell, idelta) = mean( resp_win - base_win );
    dfof_ste_noad(icell, idelta) = std( resp_win - base_win ) / sqrt(ntrial_cond);
    
end
end

figure
subplot(1,2,1)
imagesc(sig_ttest_noad); colorbar
title('visually driven by no-adapter trial')
subplot(1,2,2)
imagesc(p_ttest_noad(:,:,1)); colorbar
title('p value')
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['visual_driven_cells_noad_targ.jpg'])
close

% vis-driven by ad
sig_ttest_ad = pi * ones(ncell, ndelta); p_ttest_ad = pi * ones(ncell, ndelta);
base_avg_ad = pi * ones(ncell, ndelta);
resp_avg_ad = pi * ones(ncell, ndelta); resp_ste_ad = pi * ones(ncell, ndelta); % standard error 
dfof_avg_ad = pi * ones(ncell, ndelta); dfof_ste_ad = pi * ones(ncell, ndelta); % dF/F

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
for icell = 1 : ncell
    base_win = []; resp_win = [];
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % use only with-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 

        base_win = [base_win; mean(squeeze(tc_trial_align_ad(icell, idx, range_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trial_align_ad(icell, idx, range_resp)),2)];
    end
        
    [sig_ttest_ad(icell, idelta), p_ttest_ad(icell, idelta)] = ttest(base_win, resp_win,...
            'alpha',0.01./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg_ad(icell, idelta) = mean(base_win); % avg over trials of same ori
    resp_avg_ad(icell, idelta) = mean(resp_win);
    resp_ste_ad(icell, idelta) = std(resp_win) / sqrt(length(resp_win));

    dfof_avg_ad(icell, idelta) = mean( resp_win - base_win );
    dfof_ste_ad(icell, idelta) = std( resp_win - base_win ) / sqrt(ntrial_cond);
    
end
end

subplot(1,2,1)
imagesc(sig_ttest_ad); colorbar
title('visually driven by adapter trial')
subplot(1,2,2)
imagesc(p_ttest_ad(:,:,1)); colorbar
title('p value')
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['visual_driven_cells_adapter.jpg'])
close

save resp_noad_targ.mat dfof_avg_noad dfof_ste_noad sig_ttest_noad 
save resp_ad.mat dfof_avg_ad dfof_ste_ad sig_ttest_ad 

%% find good_fit_cells w bootstrap using new dfof

bootstrap_file = fullfile(result_folder, 'ori_across_bootstrap_runs.mat');
if ~exist(bootstrap_file)

nrun = 1000;
dfof_avg_runs = pi * ones(ncell, ndelta, nrun);
dfof_ste_runs = pi * ones(ncell, ndelta, nrun);
fit_param_runs = pi * ones(ncell, 7, nrun);
ori_pref_runs = pi * ones(ncell, nrun);

theta = deg2rad(delta_list);
disp('start bootstrap runs')
for irun = 1 : nrun
    disp(num2str(irun))

    for icell = 1 : ncell        
        for idelta = 1 : ndelta
            idx = find(delta_seq == delta_list(idelta)); 
            idx = intersect(idx, id_noad);
            
            ntrials_delta_noad = length(idx);
            bootstrap_draw = round(ntrials_delta_noad * 0.7);
            idx_run = randsample(idx, bootstrap_draw, 1); % w replacement

            % well-fit for no-ad only
            base_win = squeeze(tc_trial_align_targ(icell, idx_run, range_base));
            base_win = mean(base_win, 2); % avg over window -> [ntrial, 1]
            resp_win = squeeze(tc_trial_align_targ(icell, idx_run, range_resp));
            resp_win = mean(resp_win, 2);

            dfof_avg_runs(icell, idelta, irun) = mean( resp_win - base_win );
            dfof_ste_runs(icell, idelta, irun) = std( resp_win - base_win ) ./ sqrt(ntrials_delta_noad);
        end
        
        data = dfof_avg_runs(icell, :, irun); 
        [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
        fit_param_runs(icell, :, irun) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
    %   icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2

        ori_pref = rad2deg(u1_hat);
        ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180;
        ori_pref(ori_pref >= 180) = ori_pref(ori_pref >= 180) - 180;
        ori_pref_runs(icell, irun) = ori_pref;
    end
end
cd(result_folder)
save ori_across_bootstrap_runs.mat dfof_avg_runs dfof_ste_runs fit_param_runs ori_pref_runs

% % sanity check
% subplot(1,2,1)
% imagesc(dfof_avg_noad); colorbar
% subplot(1,2,2)
% imagesc(mean(dfof_avg_runs, 3)); colorbar
% set(gcf, 'Position', get(0, 'Screensize'));
% % imagesc(dfof_avg_runs(:,:,1))

else
    load(bootstrap_file)
end

%% get dfof & fit of [noad, 750, 250]

base_cond = cell(ncell, ndelta, ngap); resp_cond = cell(ncell, ndelta, ngap);
for icell = 1 : ncell
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % with-ad, specific isi & ori
        ntrial_cond = length(idx); 
        
        base_cond{icell, idelta, igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_base)),2); % avg over window -> [ntrial_ori, 1]
        resp_cond{icell, idelta, igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_resp)),2);
    end
end
end

sig_ttest_cond = pi * ones(ncell, ndelta, ngap); p_ttest_cond = pi * ones(ncell, ndelta, ngap);
base_avg_cond = pi * ones(ncell, ndelta, ngap);
resp_avg_cond = pi * ones(ncell, ndelta, ngap); resp_ste_cond = pi * ones(ncell, ndelta, ngap); % standard error 
dfof_avg_cond = pi * ones(ncell, ndelta, ngap); dfof_ste_cond = pi * ones(ncell, ndelta, ngap); % dF/F
for icell = 1 : ncell
for idelta = 1 : ndelta 
    for igap =  1 : ngap
        ntrial_cond = length(base_cond{icell, idelta, igap});
       [sig_ttest_cond(icell, idelta, igap), p_ttest_cond(icell, idelta, igap)] = ttest(base_cond{icell, idelta, igap}, ...
           resp_cond{icell, idelta, igap},...
           'alpha',0.01./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
            
        base_avg_cond(icell, idelta, igap) = mean(base_cond{icell, idelta, igap}); % avg over trials of same ori
        resp_avg_cond(icell, idelta, igap) = mean(resp_cond{icell, idelta, igap});
        resp_ste_cond(icell, idelta, igap) = std(resp_cond{icell, idelta, igap}) / sqrt(ntrial_cond);

        dfof_avg_cond(icell, idelta, igap) = mean( resp_cond{icell, idelta, igap} - base_cond{icell, idelta, igap} );
        dfof_ste_cond(icell, idelta, igap) = std( resp_cond{icell, idelta, igap} - base_cond{icell, idelta, igap} ) / sqrt(ntrial_cond);
    end
end
end
cd(fullfile(result_folder, 'pre-processing'))
save resp_ad_targ.mat dfof_avg_cond dfof_ste_cond

% tuning curve fit by cond
dfof_avg_750 = dfof_avg_cond(:,:,1); dfof_ste_750 = dfof_ste_cond(:,:,1); 
fit_param_750 = pi * ones(ncell, 7);
for icell = 1 : ncell
    theta = deg2rad(delta_list);
    data = dfof_avg_750(icell,:); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
    fit_param_750(icell,:) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
%   icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2
end
u1_hat_cells = fit_param_750(:,5);
ori_pref = rad2deg(u1_hat_cells);
ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180; ori_pref(ori_pref > 180) = ori_pref(ori_pref > 180) - 180;
ori_pref_cells_750 = ori_pref;

dfof_avg_250 = dfof_avg_cond(:,:,2); dfof_ste_250 = dfof_ste_cond(:,:,2);
fit_param_250 = pi * ones(ncell, 7);
for icell = 1 : ncell
    theta = deg2rad(delta_list);
    data = dfof_avg_250(icell,:); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
    fit_param_250(icell,:) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
end
u1_hat_cells = fit_param_250(:,5);
ori_pref = rad2deg(u1_hat_cells);
ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180; ori_pref(ori_pref > 180) = ori_pref(ori_pref > 180) - 180;
ori_pref_cells_250 = ori_pref;

fit_param_noad = pi * ones(ncell, 7);
for icell = 1 : ncell
    theta = deg2rad(delta_list);
    data = dfof_avg_noad(icell,:); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
    fit_param_noad(icell,:) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
end
u1_hat_cells = fit_param_noad(:,5);
ori_pref = rad2deg(u1_hat_cells);
ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180; ori_pref(ori_pref > 180) = ori_pref(ori_pref > 180) - 180;
ori_pref_cells_noad = ori_pref;

dfof_avg_merge = cat(3, dfof_avg_noad, dfof_avg_750, dfof_avg_250);
dfof_ste_merge = cat(3, dfof_ste_noad, dfof_ste_750, dfof_ste_250);
fit_param_merge = cat(3, fit_param_noad, fit_param_750, fit_param_250);
ori_pref_cells_merge = cat(2, ori_pref_cells_noad, ori_pref_cells_750, ori_pref_cells_250);

cd(result_folder)
save dfof_fit_noad_750_250.mat dfof_avg_merge dfof_ste_merge fit_param_merge ori_pref_cells_merge

%% trial trace by cond or trace noad, converted to dfof

trace_cond_dfof = cell(ncell, ndelta, ngap); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % with-ad, specific isi, specific ori
        range_trace = [1 : max(trial_len_list)]; 
        trace_cond_dfof{icell, idelta, igap} = squeeze(tc_trial_align_ad(icell, idx, range_trace)); % [ntrial, trial_len]
%         trace_cond_targ0{icell, igap} = nanmean(trace_cond_dfof{icell, 8, igap},1);
    end
end
end

vis_driven_ad = sum(sig_ttest_ad,2)>0; sum(sum(sig_ttest_ad,2)>0)
cell_list_now = find(vis_driven_ad);
trace_targ0_750 = []; trace_targ0_250 = [];

figure
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    trace_targ0_750(ii,:) = nanmean(trace_cond_dfof{icell, 8, 1},1);
    trace_targ0_250(ii,:) = nanmean(trace_cond_dfof{icell, 8, 2},1);
end
plot(mean(trace_targ0_750,1)); hold on; plot(mean(trace_targ0_250,1));
cd C:\Users\lan\Documents\repos\inter\mat
saveas(gcf, ['grand avg trace set ', num2str(iset)], 'jpg')
close

cd(result_folder)
save trace_targ0_isi.mat trace_targ0_750 trace_targ0_250

trace_noad_dfof = cell(ncell, ndelta); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    idx = intersect(id_delta, id_noad); % no-ad / merge isi, & specific ori
    range_trace = [1 : max(trial_len_list)]; 
    trace_noad_dfof{icell, idelta} = squeeze(tc_trial_align_targ(icell, idx, range_trace)); 
end
end
save trace_noad_750_250.mat trace_cond_dfof trace_noad_dfof

%% cell list by property

ncell
vis_driven_noad = sum(sig_ttest_noad,2)>0; sum(sum(sig_ttest_noad,2)>0) % ncells responsive to >= 1 targ ori: 81/97
vis_driven_ad = sum(sig_ttest_ad,2)>0; sum(sum(sig_ttest_ad,2)>0) % ncells responsive to ad: 86/97
vis_driven = vis_driven_ad | vis_driven_noad; % vis-driven cells are activated by adapter or no-ad targ
vis_driven_cell_list = find(vis_driven);
ncell_vis_driven = length(find(vis_driven))

nrun = size(ori_pref_runs,2);
ori_closeness = sort(abs(ori_pref_runs - mean(ori_pref_runs,2)), 2); % sort each row aka cell
percentile_threshold = 0.90;
percentile_idx = percentile_threshold * nrun;
ori_perc = ori_closeness(:, percentile_idx);
good_fit_cell = ori_perc < (180/ndelta);
ncell_good_fit = sum(good_fit_cell)

best_cell = vis_driven & good_fit_cell;
ncell_best = length(find(best_cell))

ori_pref_cells_noad = ori_pref_cells_merge(:,1);
pref_0_cell = find(ori_pref_cells_noad > (delta_list(end-1) + delta_list(end))/2 | ori_pref_cells_noad < delta_list(1)/2);
% pref_0_cell = pref_0_cell( dfof_avg_tg0(pref_0_cell,1)>0 & dfof_avg_tg0(pref_0_cell,2)>0); % dfof should >0
% pref_0_cell = pref_0_cell( dfof_avg_ad(pref_0_cell)>0 ); % dfof should >0
pref_0_cell = intersect(pref_0_cell, vis_driven_cell_list);

save cell_property.mat vis_driven_noad vis_driven_ad vis_driven good_fit_cell best_cell pref_0_cell

end