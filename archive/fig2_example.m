%% set up directory

close all
clear
clc

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

% need automation here. import excel
date = '200728';
mouse = '1324';
ImgFolder = '003';
time = strvcat('1308'); 
frame_rate = 30;

imouse = ['i', mouse];
result_folder = 'C:\Users\lan\Documents\repos\inter\code\200720_i1323_v3';
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' imouse];
datemouserun = [date '_' imouse '_' run_str];
fnout = fullfile(ll_fn, 'Analysis\2P\', datemouserun); 

%% load data

fName = fullfile(mworks_fn, ['data-' imouse '-' date '-' time '.mat']);
load(fName); % load behavior data "input"

CD = fullfile(data_fn, imouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata "info"

tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs_addfake.mat']); % load time course including fake targ resp
% fix bug here: tc folder in Analysis naming convention

cd C:\Users\lan\Documents\repos\inter\code\
% load ori_across_bootstrap_runs_with_replace.mat
% load ori_across_cells_cond.mat
% load resp_noad_targ.mat
% load resp_ad.mat
% load resp_ad_targ.mat
% load resp_ad_targ0deg.mat 
% load trace_by_cond.mat
% load trace_by_cond_dfof.mat

%% exp design

ntrial = input.trialSinceReset - 1; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps
% final trial discarded bc too few frames

cStimOn = double(cell2mat(input.cStimOn)); cStimOff = double(cell2mat(input.cStimOff));
cTarget = celleqel2mat_padded(input.cTargetOn); cTarget = double(cTarget);

adapter_stim_len = min(unique(cStimOff - cStimOn)) % adapter = 3-4 frames, rounded from 100 ms
isi_len = unique(cTarget - cStimOff) % ISI = 8 or 22-23 frames, rounded up from 250 or 750 ms
targ_stim_len = floor(input.targetOnTimeMs / frame_rate) % target = 3-4 frames, rounded from 100 ms
iti_len = unique(cStimOn(2:end) - (cTarget(1:end-1) + 3)) % ITI = 192-193-194 frames, 6.4-6.5 s?
trial_len_list = unique(diff(cStimOn)); % 6.9 or 7.4s

targCon = celleqel2mat_padded(input.tGratingContrast); unique(targCon); % target contrast 1
adapterCon = celleqel2mat_padded(input.tBaseGratingContrast); unique(adapterCon); % adapter contrast 0 or 1

adapterDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(adapterDir); % adapter dir === 0
ndir = length(dirs);
delta_seq = celleqel2mat_padded(input.tGratingDirectionDeg);
delta_list = unique(delta_seq); % target 8 dir (actually ori): 22.5-180. equivalent to diff from adapter
% potential bug: need to fix delta_list to convert 180 to 0
ndelta = length(delta_list); 

[nframe, ncell] = size(npSub_tc) % nframe * ncell

%% align tc by stim onset or targ onset 
% not accounting for ca latency yet. did not leave room for 1-sec-long "trial baseline"

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

%% normalize aligned tc by 1-sec-long "trial baseline"

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
% saveas(gcf, ['aligned_tc_zoomin', num2str(icell)], 'jpg')
disp('select base vs resp window accordingly, and edit below')

%% index of trials

ca_latency = 8; % monitor stim onset in frame #1 lead to neural signal in frame #9
isi_list = cTarget - cStimOff; 
ngap = length(unique(cTarget - cStimOn));

id_750 = find(isi_list > mean(isi_list)); id_750(id_750 > ntrial) = [];
id_250 = find(isi_list < mean(isi_list)); id_750(id_750 > ntrial) = [];
id_gaps = {id_750, id_250}; 
id_noad = intersect(find(targCon == 1),find(adapterCon == 0)); id_noad(id_noad > ntrial) = [];
id_ad = intersect(find(targCon == 1),find(adapterCon == 1)); id_ad(id_ad > ntrial) = [];

%% resp to no-ad targ or ad -> find vis-driven neurons

% vis-driven by noad targ
sig_ttest_noad = pi * ones(ncell, ndelta); p_ttest_noad = pi * ones(ncell, ndelta);
base_avg_noad = pi * ones(ncell, ndelta);
resp_avg_noad = pi * ones(ncell, ndelta); resp_ste_noad = pi * ones(ncell, ndelta); % standard error 
cp_win_noad = cell(ncell, ndelta);
dfof_avg_noad = pi * ones(ncell, ndelta); dfof_ste_noad = pi * ones(ncell, ndelta); % dF/F

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
for icell = 1 : ncell
    base_win = []; resp_win = [];
    
    for igap =  1 : ngap % order: 750, 250
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_noad); % use only no-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 
        
        range_adapt_base = [1:3]; 
        range_targ_resp = [9:11]; % ref section: find base & resp window in aligned tc

        base_win = [base_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_adapt_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_targ_resp)),2)];
    end
        
    [sig_ttest_noad(icell, idelta), p_ttest_noad(icell, idelta)] = ttest(base_win, resp_win,...
            'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg_noad(icell, idelta) = mean(base_win); % avg over trials of same ori
    resp_avg_noad(icell, idelta) = mean(resp_win);
    resp_ste_noad(icell, idelta) = std(resp_win) / sqrt(length(resp_win));
    cp_win_noad{icell, idelta} = [base_win, resp_win];

    dfof_avg_noad(icell, idelta) = mean( resp_win - base_win );
    dfof_ste_noad(icell, idelta) = std( resp_win - base_win ) / sqrt(ntrial_cond);
    
end
end
vis_driven_noad = sum(sig_ttest_noad,2)>0;
% % save resp_noad_targ.mat dfof_avg_noad dfof_ste_noad cp_win_noad base_avg_noad resp_avg_noad resp_ste_noad sig_ttest_noad p_ttest_noad

sum(sum(sig_ttest_noad,2)>0) % ncells responsive to >= 1 targ ori: 84/103
figure
subplot(1,2,1)
imagesc(sig_ttest_noad); colorbar
title('visually driven by no-adapter trial')
subplot(1,2,2)
imagesc(p_ttest_noad(:,:,1)); colorbar
title('p value')
set(gcf, 'Position', get(0, 'Screensize'));
% % saveas(gcf, ['visual_driven_cells_noad_targ.jpg'])
% close

%%
% vis-driven by ad
sig_ttest_ad = pi * ones(ncell, ndelta); p_ttest_ad = pi * ones(ncell, ndelta);
base_avg_ad = pi * ones(ncell, ndelta);
resp_avg_ad = pi * ones(ncell, ndelta); resp_ste_ad = pi * ones(ncell, ndelta); % standard error 
cp_win_ad = cell(ncell, ndelta);
dfof_avg_ad = pi * ones(ncell, ndelta); dfof_ste_ad = pi * ones(ncell, ndelta); % dF/F

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta));
    
for icell = 1 : ncell
    base_win = []; resp_win = [];
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % use only with-adapter trials with 1 ISI & 1 ori
        ntrial_cond = length(idx); 
        
        range_adapt_base = [1:3]; 
        range_targ_resp = [9:11]; 

        base_win = [base_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_adapt_base)),2)]; % avg over window -> [ntrial_ori, 1]
        resp_win = [resp_win; mean(squeeze(tc_trial_align_targ(icell, idx, range_targ_resp)),2)];
    end
        
    [sig_ttest_ad(icell, idelta), p_ttest_ad(icell, idelta)] = ttest(base_win, resp_win,...
            'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
    base_avg_ad(icell, idelta) = mean(base_win); % avg over trials of same ori
    resp_avg_ad(icell, idelta) = mean(resp_win);
    resp_ste_ad(icell, idelta) = std(resp_win) / sqrt(length(resp_win));
    cp_win_ad{icell, idelta} = [base_win, resp_win];

    dfof_avg_ad(icell, idelta) = mean( resp_win - base_win );
    dfof_ste_ad(icell, idelta) = std( resp_win - base_win ) / sqrt(ntrial_cond);
    
end
end
vis_driven_ad = sum(sig_ttest_ad,2)>0;
% % save resp_ad.mat dfof_avg_ad dfof_ste_ad cp_win_ad base_avg_ad resp_avg_ad resp_ste_ad sig_ttest_ad p_ttest_ad

sum(sum(sig_ttest_ad,2)>0) % ncells responsive to >= 1 targ ori: 84/103
subplot(1,2,1)
imagesc(sig_ttest_ad); colorbar
title('visually driven by no-adapter trial')
subplot(1,2,2)
imagesc(p_ttest_ad(:,:,1)); colorbar
title('p value')
set(gcf, 'Position', get(0, 'Screensize'));
% % saveas(gcf, ['visual_driven_cells_adapter.jpg'])
% close

% %% find good_fit_cells w bootstrap using new dfof
% 
% nrun = 1000;
% dfof_avg_runs = pi * ones(ncell, ndelta, nrun);
% dfof_ste_runs = pi * ones(ncell, ndelta, nrun);
% fit_param_runs = pi * ones(ncell, 7, nrun);
% ori_pref_runs = pi * ones(ncell, nrun);
% 
% theta = deg2rad(delta_list);
% disp('start bootstrap runs')
% for irun = 1 : nrun
%     disp(num2str(irun))
% 
%     for icell = 1 : ncell        
%         for idelta = 1 : ndelta
%             idx = find(delta_seq == delta_list(idelta)); 
%             idx = intersect(idx, id_noad);
%             
%             ntrials_delta_noad = length(idx);
%             bootstrap_draw = round(ntrials_delta_noad * 0.7);
%             idx_run = randsample(idx, bootstrap_draw, 1); % w replacement
% 
%             % well-fit for no-ad only
%             base_win = squeeze(tc_trial_align_targ(icell, idx_run, 1:3));
%             base_win = mean(base_win, 2); % avg over window -> [ntrial, 1]
%             resp_win = squeeze(tc_trial_align_targ(icell, idx_run, 9:11));
%             resp_win = mean(resp_win, 2);
% 
%             dfof_avg_runs(icell, idelta, irun) = mean( resp_win - base_win );
%             dfof_ste_runs(icell, idelta, irun) = std( resp_win - base_win ) ./ sqrt(ntrials_delta_noad);
%         end
%         
%         
%         data = dfof_avg_runs(icell, :, irun); 
%         [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
%         fit_param_runs(icell, :, irun) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
%     %   icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2
% 
%         ori_pref = rad2deg(u1_hat);
%         ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180;
%         ori_pref(ori_pref >= 180) = ori_pref(ori_pref >= 180) - 180;
%         ori_pref_runs(icell, irun) = ori_pref;
%     end
% end
% % % save ori_across_bootstrap_runs.mat dfof_avg_runs dfof_ste_runs fit_param_runs ori_pref_runs
% 
% % sanity check
% subplot(1,2,1)
% imagesc(dfof_avg_noad); colorbar
% subplot(1,2,2)
% imagesc(mean(dfof_avg_runs, 3)); colorbar
% set(gcf, 'Position', get(0, 'Screensize'));
% % imagesc(dfof_avg_runs(:,:,1))

%% cell list by property

load('C:\Users\lan\Documents\repos\inter\mat\V1_i1324_200728\ori_across_bootstrap_runs.mat')
ncell

vis_driven_ad = sum(sig_ttest_ad,2)>0; vis_driven_noad = sum(sig_ttest_noad,2)>0;
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

%% resp to with-adapter 0-deg targ 

sig_ttest_tg0 = pi * ones(ncell, ngap); p_ttest_tg0 = pi * ones(ncell, ngap);
base_avg_tg0 = pi * ones(ncell, ngap);
resp_avg_tg0 = pi * ones(ncell, ngap); resp_ste_tg0 = pi * ones(ncell, ngap); % standard error 
cp_win_tg0 = cell(ncell, ngap);
dfof_avg_tg0 = pi * ones(ncell, ngap); dfof_ste_tg0 = pi * ones(ncell, ngap); % dF/F

id_delta0 = find(delta_seq == 180);
for icell = 1 : ncell
    base_win = cell(1,2); resp_win = cell(1,2);
    
    for igap =  1 : ngap
        idx = intersect(intersect(id_gaps{igap}, id_delta0), id_ad); % use only with-adapter 0-deg trials
        ntrial_cond = length(idx); 
        
        range_adapt_base = [1:3]; 
        range_targ_resp = [9:11];
        base_win{1,igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_adapt_base)),2); % avg over window -> [ntrial_ori, 1]
        resp_win{1,igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_targ_resp)),2);
    end

    for igap =  1 : ngap
       [sig_ttest_tg0(icell, igap), p_ttest_tg0(icell, igap)] = ttest(base_win{1,igap}, resp_win{1,igap},...
                'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
            
        base_avg_tg0(icell, igap) = mean(base_win{1,igap}); % avg over trials of same ori
        resp_avg_tg0(icell, igap) = mean(resp_win{1,igap});
        resp_ste_tg0(icell, igap) = std(resp_win{1,igap}) / sqrt(length(resp_win{1,igap}));
        cp_win_tg0{icell, igap} = [base_win{1,igap}, resp_win{1,igap}];

        dfof_avg_tg0(icell, igap) = mean( resp_win{1,igap} - base_win{1,igap} );
        dfof_ste_tg0(icell, igap) = std( resp_win{1,igap} - base_win{1,igap} ) / sqrt(ntrial_cond);
    end
end
% % save resp_ad_targ0.mat dfof_avg_tg0 dfof_ste_tg0 cp_win_tg0 base_avg_tg0 resp_avg_tg0 resp_ste_tg0 sig_ttest_tg0 p_ttest_tg0

%% resp to with-ad targ by cond (dir & isi)

base_cond = cell(ncell, ndelta, ngap); resp_cond = cell(ncell, ndelta, ngap);
for icell = 1 : ncell
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % with-ad, specific isi & ori
        ntrial_cond = length(idx); 
        
        range_adapt_base = [1:3]; 
        range_targ_resp = [9:11];
        base_cond{icell, idelta, igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_adapt_base)),2); % avg over window -> [ntrial_ori, 1]
        resp_cond{icell, idelta, igap} = mean(squeeze(tc_trial_align_targ(icell, idx, range_targ_resp)),2);
    end
end
end

sig_ttest_cond = pi * ones(ncell, ndelta, ngap); p_ttest_cond = pi * ones(ncell, ndelta, ngap);
base_avg_cond = pi * ones(ncell, ndelta, ngap);
resp_avg_cond = pi * ones(ncell, ndelta, ngap); resp_ste_cond = pi * ones(ncell, ndelta, ngap); % standard error 
cp_win_cond = cell(ncell, ndelta, ngap);
dfof_avg_cond = pi * ones(ncell, ndelta, ngap); dfof_ste_cond = pi * ones(ncell, ndelta, ngap); % dF/F
for icell = 1 : ncell
for idelta = 1 : ndelta 
    for igap =  1 : ngap
        ntrial_cond = length(base_cond{icell, idelta, igap});
       [sig_ttest_cond(icell, idelta, igap), p_ttest_cond(icell, idelta, igap)] = ttest(base_cond{icell, idelta, igap}, ...
           resp_cond{icell, idelta, igap},...
           'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
            
        base_avg_cond(icell, idelta, igap) = mean(base_cond{icell, idelta, igap}); % avg over trials of same ori
        resp_avg_cond(icell, idelta, igap) = mean(resp_cond{icell, idelta, igap});
        resp_ste_cond(icell, idelta, igap) = std(resp_cond{icell, idelta, igap}) / sqrt(ntrial_cond);
        cp_win_cond{icell, idelta, igap} = [base_cond{icell, idelta, igap}, resp_cond{icell, idelta, igap}];

        dfof_avg_cond(icell, idelta, igap) = mean( resp_cond{icell, idelta, igap} - base_cond{icell, idelta, igap} );
        dfof_ste_cond(icell, idelta, igap) = std( resp_cond{icell, idelta, igap} - base_cond{icell, idelta, igap} ) / sqrt(ntrial_cond);
    end
end
end
% % save resp_ad_targ.mat dfof_avg_cond dfof_ste_cond cp_win_cond base_avg_cond resp_avg_cond resp_ste_cond sig_ttest_cond p_ttest_cond

%% tuning curve fit by cond

dfof_avg_750 = dfof_avg_cond(:,:,1);
dfof_ste_750 = dfof_ste_cond(:,:,1);
cp_win_750 = cp_win_cond(:,:,1);
sig_ttest_750 = sig_ttest_cond(:,:,1);
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

dfof_avg_250 = dfof_avg_cond(:,:,2);
dfof_ste_250 = dfof_ste_cond(:,:,2);
cp_win_250 = cp_win_cond(:,:,2);
sig_ttest_250 = sig_ttest_cond(:,:,2);
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

% % save ori_across_cells_cond.mat dfof_avg_merge dfof_ste_merge fit_param_merge ori_pref_cells_merge

%% Fig 1E: compare ad_resp vs 0-deg targ_resp_750/250 of vis_driven cells
% buggy: why is resp_tg/resp_ad larger than 1???

ori_pref_cells_noad = ori_pref_cells_merge(:,1);
pref_0_cell = find(ori_pref_cells_noad > (delta_list(end-1) + delta_list(end))/2 | ori_pref_cells_noad < delta_list(1)/2);
pref_0_cell = pref_0_cell( dfof_avg_tg0(pref_0_cell,1)>0 & dfof_avg_tg0(pref_0_cell,2)>0); % dfof should >0
pref_0_cell = pref_0_cell( dfof_avg_ad(pref_0_cell)>0 ); % dfof should >0
pref_0_cell = intersect(pref_0_cell, vis_driven_cell_list);

[count, idx] = histc(dfof_avg_ad(pref_0_cell), 0:0.05:ceil(max(dfof_avg_ad(pref_0_cell))*10)/10);
count_nonzero_id = find(count~=0);
% histogram(dfof_avg_ad(pref_0_cell), 8)
resp_bin_ad = accumarray(idx(:), dfof_avg_ad(pref_0_cell),[],@mean)
resp_bin_tg = [];
resp_bin_std = []; 
resp_bin_ste = zeros(length(count),2);
resp_bin_tg(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,1),[],@mean);
resp_bin_tg(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,2),[],@mean)
resp_bin_std(:,1) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,1),[],@std);
resp_bin_std(:,2) = accumarray(idx(:),dfof_avg_tg0(pref_0_cell,2),[],@std);
% resp_bin_ste(count~=0, :) = resp_bin_std(count~=0, :) ./ sqrt(count(count~=0)); % use std or ste?

figure
scatter(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, 'b.'); hold on 
scatter(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, 'r.')
errorbar(resp_bin_ad, resp_bin_tg(:,1)./resp_bin_ad, resp_bin_std(:,1), 'b', 'LineStyle','none');
errorbar(resp_bin_ad, resp_bin_tg(:,2)./resp_bin_ad, resp_bin_std(:,2), 'r', 'LineStyle','none');
for itext = 1 : length(count_nonzero_id)
    text(resp_bin_ad(count_nonzero_id(itext)), ...
        resp_bin_tg(count_nonzero_id(itext),2)./resp_bin_ad(count_nonzero_id(itext)) - resp_bin_std(count_nonzero_id(itext),2) - 0.1, ...
        ['n=', num2str(count(count_nonzero_id(itext)))], 'HorizontalAlignment', 'center')
end
ylim([0,3])
legend('isi 750', 'isi 250')
% % saveas(gcf, ['Fig 1E cells prefer 0 deg binned by ad-resp.jpg'])
% close

%% Fig 1C (approx): time course of adaptation recovery for all cells
% switched to median bc of super negative outlier

isi_sec = [4, 0.750, 0.250];
norm_targ_resp_mean = [1, median((dfof_avg_tg0(vis_driven,1)./dfof_avg_ad(vis_driven))), ...
    median(dfof_avg_tg0(vis_driven,2)./dfof_avg_ad(vis_driven))];

data = table(isi_sec', norm_targ_resp_mean');
modelfun = @(param, data) 1 - param(1) * exp(param(2) * data(:, 1)); % Y = 1 - y0 * exp(k*x)
param_init = [0.75, -2]; % init guess values 
mdl = fitnlm(data, modelfun, param_init);
param = mdl.Coefficients{:, 'Estimate'}
x = 0:0.01:4
yfit_decay = 1 - param(1) * exp(param(2) * x);

scatter(isi_sec, norm_targ_resp_mean)
hold on;
plot(x, yfit_decay, 'r-');
grid on; grid minor
% formulaString = sprintf('Y = 1 - %.3f * exp(%.3f * X)', param_decay(1), param_decay(2))
xlim([0, 4.5])
ylim([0, 1.3])
xlabel('ISI (s)')
ylabel('normalized dF/F')
legend off
% % saveas(gcf, ['Fig 1C time course of adaptation recovery of vis-driven cells.jpg'])
% close

%% trial trace by cond or trace noad, converted to dfof

trace_cond_dfof = cell(ncell, ndelta, ngap); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    
    for igap =  1 : ngap 
        idx = intersect(intersect(id_gaps{igap}, id_delta), id_ad); % with-ad, specific isi, specific ori
        range_trace = [1 : max(trial_len_list)]; 
        trace_cond_dfof{icell, idelta, igap} = squeeze(tc_trial_align_ad(icell, idx, range_trace)); % [ntrial, trial_len]
    end
end
end

trace_noad_dfof = cell(ncell, ndelta); 
for icell = 1 : ncell    
for idelta = 1 : ndelta 
    id_delta = find(delta_seq == delta_list(idelta));
    idx = intersect(id_delta, id_noad); % no-ad / merge isi, & specific ori
    range_trace = [1 : max(trial_len_list)]; 
    trace_noad_dfof{icell, idelta} = squeeze(tc_trial_align_targ(icell, idx, range_trace)); 
end
end

% % save trace_by_cond_dfof.mat trace_cond_dfof trace_noad_dfof

%% Fig 2B: cell resp by dir & condition (no_ad, 750, 250) for "best cells"

cell_list_now = vis_driven_cell_list;
for ii = 27 % 1 : length(cell_list_now)
icell = cell_list_now(ii);
trace_noad_avg = []; trace_cond_avg_750 = []; trace_cond_avg_250 = [];

for idelta = 1 : ndelta
    trace_noad_avg(idelta, :) = nanmean(trace_noad_dfof{icell, idelta}, 1);
    trace_cond_avg_750(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 1}, 1);
    trace_cond_avg_250(idelta, :) = nanmean(trace_cond_dfof{icell, idelta, 2}, 1);
end

adjusted_trace_len = 3.5 * frame_rate; % trace len 3.5 s according to Jin2019 Fig 2B
trace_noad_avg = [trace_noad_avg(end,:); trace_noad_avg]; % wrap around
trace_noad_avg = trace_noad_avg(:, 1:adjusted_trace_len); % trim
trace_cond_avg_750 = [trace_cond_avg_750(end,:); trace_cond_avg_750]; 
trace_cond_avg_750 = trace_cond_avg_750(:, 1:adjusted_trace_len);
trace_cond_avg_250 = [trace_cond_avg_250(end,:); trace_cond_avg_250]; 
trace_cond_avg_250 = trace_cond_avg_250(:, 1:adjusted_trace_len);

x = [length(trace_noad_avg), length(trace_cond_avg_750), length(trace_cond_avg_250)];
xmax = max(x);
y = [trace_noad_avg, trace_cond_avg_750, trace_cond_avg_250];
ymin = min(y(:));
ymax = max(y(:)) + 0.05 * (max(y(:)) - min(y(:)));

figure('units','normalized','outerposition',[0 0 1/2 1]);
for col = 1 : 3
    for idelta = 1 : (ndelta+1)
        subplot(ndelta+1, 3, col + 3*(idelta-1));
        if col == 1, plot(trace_noad_avg(idelta, :)) 
            if idelta == 1, title('no adapter'), end
        elseif col == 2, plot(trace_cond_avg_750(idelta, :))
             if idelta == 1, title('isi 750'), end
        elseif col == 3, plot(trace_cond_avg_250(idelta, :))
             if idelta == 1, title('isi 250'), end
        end
%         axis off
        xlim([0, xmax])
        ylim([ymin, ymax])
        xticks(0 : 30 : xmax)
        yticks(round(ymin*10)/10 : round((ymax-ymin)*10)/20 : round(ymax*10)/10)
    end
end
% % saveas(gcf, ['dfof trace ', num2str(icell), '.jpg'])
% close    
end

%% adp san check
% suspicious: dfof_equiv_noad avg get super close to 0, bc it contains lots of negative dfof
% is this normal???

% with-adapter / no-adapter resp to same targ with same isi
cell_list_now = pref_0_cell;
adp_ratio = zeros(length(cell_list_now), ndelta, ngap);

for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
% for icell = 1:ncell
    for idelta = 1:ndelta
        id_delta = find(delta_seq == delta_list(idelta));
        for igap = 1:ngap
            idx_now_ad = intersect(intersect(id_gaps{igap}, id_delta), id_ad);
            idx_now_noad = intersect(id_delta, id_noad);
            dfof_equiv_ad = mean(squeeze(tc_trial_align_targ(icell, idx_now_ad, 9:11)),2) - mean(squeeze(tc_trial_align_targ(icell, idx_now_ad, 1:3)),2);
            dfof_equiv_noad = mean(squeeze(tc_trial_align_targ(icell, idx_now_noad, 9:11)),2) - mean(squeeze(tc_trial_align_targ(icell, idx_now_noad, 1:3)),2);
            adp_ratio(ii, idelta, igap) = mean(dfof_equiv_ad) / mean(dfof_equiv_noad) - 1;
%             adp_ratio(ii, idelta, igap) = mean(dfof_equiv_ad(dfof_equiv_ad>0)) / mean(dfof_equiv_noad(dfof_equiv_noad>0)) - 1;
        end
    end
end

all = sum(adp_ratio(:)>-Inf)
outlier_facil = find(adp_ratio > mean(adp_ratio(:) + 3*std(adp_ratio(:))) | adp_ratio < mean(adp_ratio(:) - 3*std(adp_ratio(:))))
facilitated = sum(adp_ratio(:)>0)
inhibited = sum(adp_ratio(:)<=0)

adp_ratio(outlier_facil) = NaN;
nanmean(adp_ratio(:))
nanstd(adp_ratio(:))

t750 = squeeze(adp_ratio(:,:,1));
t250 = squeeze(adp_ratio(:,:,2));
errorbar(nanmean(t750, 1), nanstd(t750, 1), 'b'); hold on;
errorbar(nanmean(t250, 1), nanstd(t250, 1), 'r');
line([0,9], [0, 0], 'Color', 'g', 'LineWidth', 1);
xlabel('ori');ylabel('adp ratio')

figure
subplot(1,2,1)
histogram(adp_ratio(:), 30)
yl = ylim;
line([0, 0], [0, yl(2)], 'Color', 'g', 'LineWidth',1, 'LineStyle','--');
xlabel('adp index')
ylabel('count across cell ori isi')


% check only with-ad targ0 vs ad || with-ad targ0 vs no-ad targ0
cell_list_now = pref_0_cell;
adp_ratio_targ0 = zeros(length(cell_list_now), ngap);

for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
% for icell = 1:ncell
    for idelta = 8 % targ0 only! adp is ori-specific
        id_delta = find(delta_seq == delta_list(idelta));
        for igap = 1:ngap
            idx_now_ad_targ = intersect(intersect(id_gaps{igap}, id_delta), id_ad);
%             idx_now_noad_targ = intersect(id_delta, id_noad);
            dfof_equiv_ad_targ = mean(squeeze(tc_trial_align_targ(icell, idx_now_ad_targ, 9:11)),2)...
                               - mean(squeeze(tc_trial_align_targ(icell, idx_now_ad_targ, 1:3)),2);
            dfof_equiv_ad = mean(squeeze(tc_trial_align_ad(icell, idx_now_ad_targ, 9:11)),2)...
                          - mean(squeeze(tc_trial_align_ad(icell, idx_now_ad_targ, 1:3)),2);
%             dfof_equiv_noad_targ = mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, 9:11)),2) - mean(squeeze(tc_trial_align_targ(icell, idx_now_noad_targ, 1:3)),2);
%             adp_ratio_targ0(ii, igap) = mean(dfof_equiv_ad_targ(dfof_equiv_ad_targ>0)) / mean(dfof_equiv_ad(dfof_equiv_ad>0)) - 1;
            adp_ratio_targ0(ii, igap) = mean(dfof_equiv_ad_targ) / mean(dfof_equiv_ad) - 1;
        end
    end
end

all = sum(adp_ratio_targ0(:)>-Inf)
outlier_facil0 = find(adp_ratio_targ0 > mean(adp_ratio_targ0(:) + 3*std(adp_ratio_targ0(:))))
facilitated = sum(adp_ratio_targ0(:)>0)
inhibited = sum(adp_ratio_targ0(:)<=0)

adp_ratio_targ0(outlier_facil0) = NaN;
nanmean(adp_ratio_targ0, 1)
nanstd(adp_ratio_targ0, 1)

subplot(1,2,2)
histogram(adp_ratio_targ0(:), 10)
yl = ylim;
line([0, 0], [0, yl(2)], 'Color', 'g', 'LineWidth',1, 'LineStyle','--');
xlabel('adp index')
ylabel('count across cell ori isi')

% %% debug 2C
% 
% t_noad = mean(dfof_avg_merge(:,:,1),1);
% t_750 = mean(dfof_avg_merge(:,:,2),1);
% t_250 = mean(dfof_avg_merge(:,:,3),1);
% plot(t_noad, 'r'); hold on
% plot(t_750, 'g');
% plot(t_250, 'b');
% legend('noad', '750', '250', 'Location', 'southwest')

%% Fig 2C: tuning curve fit by condition (noad / 750 / 250) for vis-driven cells

theta_finer = deg2rad(0:1:179);
subplot_title = {'control', 'isi 750 ms', 'isi 250 ms'};

% cell_list_now = vis_driven_cell_list;
cell_list_now = find(good_fit_cell);
for ii = 1 %: length(cell_list_now)
    icell = cell_list_now(ii);
    figure('units','normalized','outerposition',[0 0 1/2 1]);
    
for row = 1 : 3
    subplot(3,1,row)
    
    dfof_avg_now = dfof_avg_merge(:,:,row);
    dfof_ste_now = dfof_ste_merge(:,:,row);
    fit_param_now = fit_param_merge(:,:,row);
    
    errorbar([0,delta_list], [dfof_avg_now(icell,end), dfof_avg_now(icell,:)], ...
        [dfof_ste_now(icell,end), dfof_ste_now(icell,:)], 'LineStyle','none')
    hold on
    scatter([0,delta_list], [dfof_avg_now(icell,end), dfof_avg_now(icell,:)], 'b')

    t = num2cell(fit_param_now(icell, 2:end)); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = deal(t{:});    
    y_fit(row,:) = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(theta_finer - u1_hat))-1));
    plot(rad2deg(theta_finer), y_fit(row,:), 'LineWidth', 1)
    
    ori_pref = ori_pref_cells_merge(icell, row);
    scatter(ori_pref, b_hat + R1_hat, 'r*') % mark pref ori of fit
    
    xlim([0-10, 180+10])
    tempmin = dfof_avg_merge - dfof_ste_merge; tempmin = tempmin(icell, :, :); ymin = min(tempmin(:));
    tempmax = dfof_avg_merge + dfof_ste_merge; tempmax = tempmax(icell, :, :); ymax = max(tempmax(:));
    padding = (ymax - ymin) ./ 50;
    ylim([ymin - padding, ymax + padding])
    
    xlabel('orientation (deg)')
    ylabel('dF/F')
    title(subplot_title{row})
end
%     % saveas(gcf, ['ori tuning across cond cell ', num2str(icell)], 'jpg')
%     close
end

%% Fig 2D: targ degree distance from adapter changes dfof resp
% used median instead of avg, bc some dfof_dis_noad are too close to 0, and dragged dfof_dis_norm super high
% bug: how to get rid of large outlier when using avg?

dis_seq = delta_seq; 
dis_seq(dis_seq > 90) = 180 - dis_seq(dis_seq > 90); 
dis_list = unique(dis_seq);

% mark trial number of each delta
[C,ia,ic] = unique(delta_seq);
for idelta = 1 : length(delta_list)
    ntrial_delta(idelta) = length(find(ic == idelta));
end
dis_deltas = {[8], [1,7], [2,6], [3,5], [4]};

dfof_dis_noad = {}; dfof_dis_targ = {};
cell_list_now = find(vis_driven_noad); 
% cell_list_now = find(good_fit_cell);
ncell_now = length(cell_list_now);
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    dfof_avg_noad_now = squeeze(dfof_avg_merge(icell, :, 1));
    dfof_avg_noad_targ_max = max(dfof_avg_merge(icell, :, 1)); 
    dfof_avg_noad_targ0 = squeeze(dfof_avg_merge(icell, 8, 1));  % use only noad targ 0 to normlz

for igap = 1 : ngap
    dfof_avg_isi_now = squeeze(dfof_avg_merge(icell, :, igap+1)); % 750&250

    for idis = 1 : length(dis_list)
        id_dis_delta = dis_deltas{idis};
        
        dfof_dis_noad{ii, idis, igap} = dfof_avg_noad_now(id_dis_delta) .* ntrial_delta(id_dis_delta);
        dfof_dis_noad{ii, idis, igap}= sum(dfof_dis_noad{ii, idis, igap}) ./ sum(ntrial_delta(id_dis_delta));
%         dfof_dis_noad{ii, idis, igap} = mean(); % sanity check using no-weight avg
        dfof_dis_targ{ii, idis, igap}= dfof_avg_isi_now(id_dis_delta) .* ntrial_delta(id_dis_delta);
        dfof_dis_targ{ii, idis, igap}= sum(dfof_dis_targ{ii, idis, igap}) ./ sum(ntrial_delta(id_dis_delta));
        
%         dfof_dis_norm(ii, idis, igap) = (dfof_dis_targ{ii, idis, igap} - dfof_dis_noad{ii, idis, igap})...
%            ./ dfof_dis_noad{ii, idis, igap};
        dfof_dis_norm(ii, idis, igap) = (dfof_dis_targ{ii, idis, igap} - dfof_dis_noad{ii, idis, igap})...
           ./ dfof_avg_noad_targ_max;
    end
end
end

for igap = 1 : ngap
    for idis = 1 : length(dis_list)
        dfof_dis_norm_avg(idis, igap) = mean(dfof_dis_norm(:, idis, igap));
        dfof_dis_norm_median(idis, igap) = median(dfof_dis_norm(:, idis, igap));
        dfof_dis_norm_ste(idis, igap) = std(dfof_dis_norm(:, idis, igap))./size(dfof_dis_norm, 1);
    end
end

color_list = {[0,0,1], [1,0,0]};
figure
for igap = 1 : ngap
    hold on
    scatter(dis_list, dfof_dis_norm_avg(:, igap))
%     scatter(dis_list, dfof_dis_norm_median(:, igap))
end
for igap = 1 : ngap
    hold on
    errorbar(dis_list, dfof_dis_norm_avg(:, igap), dfof_dis_norm_ste(:, igap), 'color', color_list{igap}) %, 'LineStyle','none')
%     errorbar(dis_list, dfof_dis_norm_median(:, igap), dfof_dis_norm_ste(:, igap), 'color', color_list{igap})
end
line([0-5, 180+5], [0, 0], 'Color', 'g', 'LineWidth', 1);
% ylim([-1, 0.2])
xlim([0-5, 90+5])
legend('isi 750', 'isi 250', 'Location', 'southeast'); legend boxoff
xlabel('|Test - Adapter| (deg)')
ylabel('delta norm dF/F')
title(['n = ', num2str(ncell_now), ' cells'])
% % saveas(gcf, ['Fig 2D response changes w targ-ad distance.jpg'])
% close

% % mark trial number of each distance
% [C,ia,ic] = unique(dis_seq);
% for idis = 1 : length(dis_list)
%     ntrial_dis(idis) = length(find(ic == idis));
% end
% for itext = 1 : length(dis_list)
%     text(dis_list(itext), ...
%        dfof_dis_norm_avg(itext, 1) + dfof_dis_norm_ste(itext, 1) + 0.02, ...
%         ['n=', num2str(ntrial_dis(itext))], 'HorizontalAlignment', 'center')
% end

%% Fig 2E: resp changes w |ori_pref - ori_ad| distance
% very buggy
% restrict to well-fitted cell? 

ori_pref_binned = ori_pref_cells_merge(:,1); % pref measured by noad control
ori_pref_binned(ori_pref_binned > 90) = 180 - ori_pref_binned(ori_pref_binned > 90); 
ori_pref_binned(ori_pref_binned<20) = 0;
ori_pref_binned(ori_pref_binned>=20 & ori_pref_binned<=70) = 45;
ori_pref_binned(ori_pref_binned>70) = 90;
% figure; histogram(ori_pref_binned)

cell_list_now = find(vis_driven & good_fit_cell);
dfof_peak_norm = []; dfof_noad_peak = []; dfof_isi_peak = [];
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    
    t = num2cell(fit_param_merge(icell, 2:end, 1)); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = deal(t{:});
    dfof_noad_peak(ii) = b_hat + R1_hat;

for igap = 1 : ngap
    t = num2cell(fit_param_merge(icell, 2:end, igap+1)); % 750,250
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = deal(t{:});
    dfof_isi_peak(ii, igap) = b_hat + R1_hat;

    dfof_peak_norm(ii, igap) = (dfof_isi_peak(ii, igap))./ dfof_noad_peak(ii);
end
end

ori_pref_binned_list = unique(ori_pref_binned);
dfof_peak_norm_avg = []; dfof_peak_norm_ste = []; ncell_dis = [];
for igap = 1 : ngap
    for idis = 1 : length(ori_pref_binned_list)
        dis_now = ori_pref_binned(vis_driven & good_fit_cell) == ori_pref_binned_list(idis);
        ncell_dis(idis) = sum(dis_now);
        dfof_peak_norm_avg(idis, igap) = mean(dfof_peak_norm(dis_now, igap));
        dfof_peak_norm_ste(idis, igap) = std(dfof_peak_norm(dis_now, igap)) ./ length(dfof_peak_norm(dis_now, igap));
    end
end

color_list = {[0,0,1], [1,0,0]};
figure
for igap = 1 : ngap
    hold on
    scatter(ori_pref_binned_list, dfof_peak_norm_avg(:, igap))
end
for igap = 1 : ngap
    hold on
    errorbar(ori_pref_binned_list, dfof_peak_norm_avg(:, igap), dfof_peak_norm_ste(:, igap),...
        'color', color_list{igap}) %, 'LineStyle','none')
end
for itext = 1 : length(ori_pref_binned_list)
    text(ori_pref_binned_list(itext), ...
       dfof_peak_norm_avg(itext, 1) + dfof_peak_norm_ste(itext, 1) + 0.05, ...
        ['n=', num2str(ncell_dis(itext))], 'HorizontalAlignment', 'center')
end

line([0-5, 180+5], [1, 1], 'Color', 'g', 'LineWidth', 1);
xlim([0-5, 90+5])
% ylim([0, 1.1])
legend('isi 750', 'isi 250', 'Location','southeast')
xlabel('|Pref - Adapter| (deg)')
ylabel('norm pred peak resp')
grid on; grid minor
% % saveas(gcf, ['Fig 2E response changes w pref-ad distance.jpg'])
% close

%% Fig 2F: ori_pref changes w |ori_pref - ori_ad| distance
% had to use median instead of mean so that plot makes sense...

ori_pref_dist = ori_pref_cells_merge;
ori_pref_dist(ori_pref_dist > 90) = 180 - ori_pref_dist(ori_pref_dist > 90); 

cell_list_now = find(vis_driven & good_fit_cell);
ori_pref_shift = pi *ones(length(cell_list_now), ngap);
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    dfof_noad_pref(ii) = ori_pref_dist(ii, 1);

for igap = 1 : ngap
    dfof_isi_pref(ii, igap) = ori_pref_dist(ii, igap+1);
    ori_pref_shift(ii, igap) = dfof_isi_pref(ii, igap) - dfof_noad_pref(ii);
end
end


ori_pref_shift_avg = []; ori_pref_shift_ste = []; ncell_dis = [];
for igap = 1 : ngap
    for idis = 1 : length(ori_pref_binned_list)
        dis_now = ori_pref_binned(vis_driven & good_fit_cell) == ori_pref_binned_list(idis);
        ncell_dis(idis) = sum(dis_now);
        ori_pref_shift_avg(idis, igap) = median(ori_pref_shift(dis_now, igap));
        ori_pref_shift_ste(idis, igap) = std(ori_pref_shift(dis_now, igap)) ./ length(ori_pref_shift(dis_now, igap));
    end
end

color_list = {[0,0,1], [1,0,0]};
figure; hold on
for igap = 1 : ngap    
    scatter(ori_pref_binned_list, ori_pref_shift_avg(:, igap))
end
for igap = 1 : ngap
    errorbar(ori_pref_binned_list, ori_pref_shift_avg(:, igap), ori_pref_shift_ste(:, igap),...
        'color', color_list{igap}) %, 'LineStyle','none')
end
for itext = 1 : length(ori_pref_binned_list)
    text(ori_pref_binned_list(itext), ...
       ori_pref_shift_avg(itext, 2) + ori_pref_shift_ste(itext, 2) + 1, ...
        ['n=', num2str(ncell_dis(itext))], 'HorizontalAlignment', 'center')
end

line([0-5, 180+5], [0, 0], 'Color', 'g', 'LineWidth', 1);
xlim([0-5, 90+5])
ylim([-10, 20])
legend('isi 750', 'isi 250', 'Location','southeast')
% % saveas(gcf, ['Fig 2F pref ori changes w pref-ad distance.jpg'])
% close

%% Fig 2G: OSI changes w |ori_pref - ori_ad| distance
% why osi decrease after adp? poor fit for 750/250?

cell_list_now = find(vis_driven & good_fit_cell);
osi_shift = [];
for ii = 1 : length(cell_list_now)
    icell = cell_list_now(ii);
    
    t = num2cell(fit_param_merge(icell, 2:end, 1)); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = deal(t{:});
%     y_fit(row,:) = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(theta_finer - u1_hat))-1));
    resp_pref = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(u1_hat - u1_hat))-1));
    resp_orth = max(0, b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(u1_hat + pi/2 - u1_hat))-1)));
    osi_noad(ii) = (resp_pref - resp_orth) / (resp_pref + resp_orth);

for igap = 1 : ngap
    
    t = num2cell(fit_param_merge(icell, 2:end, 1+igap)); 
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = deal(t{:});
    resp_pref = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(u1_hat - u1_hat))-1));
    resp_orth = max(0, b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(u1_hat + pi/2 - u1_hat))-1)));
    
    osi_isi(ii, igap) = (resp_pref - resp_orth) / (resp_pref + resp_orth);
    osi_shift(ii, igap) = osi_isi(ii, igap) - osi_noad(ii);
end
end


osi_shift_avg = []; osi_shift_ste = []; ncell_dis = [];
for igap = 1 : ngap
    for idis = 1 : length(ori_pref_binned_list)
        dis_now = ori_pref_binned(vis_driven & good_fit_cell) == ori_pref_binned_list(idis);
        ncell_dis(idis) = sum(dis_now);
        osi_shift_avg(idis, igap) = mean(osi_shift(dis_now, igap));
        osi_shift_ste(idis, igap) = std(osi_shift(dis_now, igap)) ./ length(osi_shift(dis_now, igap));
    end
end

color_list = {[0,0,1], [1,0,0]};
figure; hold on
for igap = 1 : ngap    
    scatter(ori_pref_binned_list, osi_shift_avg(:, igap))
end
for igap = 1 : ngap
    errorbar(ori_pref_binned_list, osi_shift_avg(:, igap), osi_shift_ste(:, igap),...
        'color', color_list{igap}) %, 'LineStyle','none')
end
for itext = 1 : length(ori_pref_binned_list)
    text(ori_pref_binned_list(itext), ...
       osi_shift_avg(itext, 1) - osi_shift_ste(itext, 1) - 0.01, ...
        ['n=', num2str(ncell_dis(itext))], 'HorizontalAlignment', 'center')
end

xlim([0-5, 90+5])
legend('isi 750', 'isi 250', 'Location','southeast')
line([0-5, 180+5], [0, 0], 'Color', 'g', 'LineWidth', 1);
