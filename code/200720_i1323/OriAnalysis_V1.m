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
date = '200720';
ImgFolder = '003';
time = strvcat('1133'); 
mouse = 'i1323';
frame_rate = 30;

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
fnout = fullfile(ll_fn, 'Analysis\2P\', datemouserun); 

%% load data

fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName); % load behavior data "input"

CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata "info"

mouse = '1323';
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs_addfake.mat']); % load time course including fake targ resp

%% exp design

ntrial = input.trialSinceReset; % 464 = 8 dir * 2 adapter contrast * 2 ISI * 14.5 reps

cStimOn = cell2mat(input.cStimOn);
cStimOff = cell2mat(input.cStimOff);
cTarget = celleqel2mat_padded(input.cTargetOn); cTarget = int64(cTarget);

adapter_stim_len = unique(cStimOff - cStimOn) % adapter = 3-4 frames, rounded from 100 ms
isi_len = unique(cTarget - cStimOff) % ISI = 8 or 22-23 frames, rounded up from 250 or 750 ms
targ_stim_len = floor(input.targetOnTimeMs / frame_rate) % target = 3-4 frames, rounded from 100 ms
iti_len = unique(cStimOn(2:end) - (cTarget(1:end-1) + 3)) % ITI = 192-193-194 frames, 6.4-6.5 s?

targCon = celleqel2mat_padded(input.tGratingContrast);
unique(targCon); % target contrast 1
adapterCon = celleqel2mat_padded(input.tBaseGratingContrast);
unique(adapterCon); % adapter contrast 0 or 1
% ind_noadapter = intersect(find(targCon == 1),find(adapterCon == 0));

adapterDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(adapterDir); % adapter dir === 0
ndir = length(dirs);
delta_seq = celleqel2mat_padded(input.tGratingDirectionDeg);
delta_list = unique(delta_seq); % target 8 dir (actually ori): 22.5-180. equivalent to diff from adapter
ndelta = length(delta_list); 

[nframe, ncell] = size(npSub_tc) % nframe * ncell

%% align by trial start

trial_len = diff(cStimOn);
trial_len_final = nframe - cStimOn(end);
trial_len = [trial_len, trial_len_final];

tc_trials = zeros(ncell, ntrial, max(trial_len));
for icell = 1:ncell
    npSub_tc_cell = npSub_tc(:,icell);
    for itrial = 1:ntrial
        start_id = cStimOn(itrial);
        tc_trials(icell, itrial, :) = [npSub_tc_cell(start_id : start_id + trial_len(itrial) - 1); ...
            NaN(max(trial_len) - trial_len(itrial), 1)];
    end
end

% % sanity check
% size(tc_trials) % ncell * ntrial * trial_len (padded by nan)
% temp = tc_trials(:);
% sum(isnan(temp))/length(temp)
% 
% trial_len = double(trial_len);
% unique_len = unique(trial_len)
% n = histc(trial_len, unique_len);
% n_nan = sum((max(unique_len) - unique_len) .* n) .* ncell
% sum(isnan(temp))

%% cells sensitive to target dir
ca_latency = 7; % Ca signal latency: signal from frame #1 shows up at frame #8
target_relative = cTarget - cStimOn; % unique(cTarget - cStimOn) = [11 26]
targ_start = 1 + target_relative + ca_latency; % time to receive targ resp signal

sig_ttest = pi * ones(ncell, ndelta, ngap); p_ttest = pi * ones(ncell, ndelta, ngap);
base_avg = pi * ones(ncell, ndelta, ngap);
resp_avg = pi * ones(ncell, ndelta, ngap); resp_ste = pi * ones(ncell, ndelta, ngap); % standard error 
dfof_avg = pi * ones(ncell, ndelta, ngap); dfof_ste = pi * ones(ncell, ndelta, ngap); % dF/F

targ_start_list = unique(targ_start);
ngap = length(targ_start_list);
for igap =  1 : ngap % ntrial per isi is equal but not distributed evenly to every delta
    id_targ = find(targ_start == targ_start_list(igap));
    range_adapt_base = [targ_start_list(igap) - targ_stim_len : targ_start_list(igap) - 1]; % adapted baseline just bef targ onset
    range_targ_resp = [targ_start_list(igap) : targ_start_list(igap) + targ_stim_len - 1]; % targ onset til targ fin

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta)); 
    idx = intersect(id_targ, id_delta);
    ntrial_cond = length(idx); % condition = specific isi (targ onset frame#) & targ dir
    
    for icell = 1 : ncell

        base_win = squeeze(tc_trials(icell, idx, range_adapt_base));
        base_win = mean(base_win, 2); % avg over window -> [ntrial_ori, 1]
        resp_win = squeeze(tc_trials(icell, idx, range_targ_resp));
        resp_win = mean(resp_win, 2);

        [sig_ttest(icell, idelta, igap), p_ttest(icell, idelta, igap)] = ttest(base_win, resp_win,...
            'alpha',0.05./(ntrial_cond - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
        base_avg(icell, idelta, igap) = mean(base_win); % avg over trials of same ori
        resp_avg(icell, idelta, igap) = mean(resp_win);
        resp_ste(icell, idelta, igap) = std(resp_win) / sqrt(length(resp_win));

        dfof_avg(icell, idelta, igap) = mean( (resp_win - base_win) ./ mean(base_win) );
        dfof_ste(icell, idelta, igap) = std( (resp_win - base_win) ./ mean(base_win) ) / sqrt(ntrial_cond);
    end
end
end

t = sig_ttest(:,:,1) + sig_ttest(:,:,2);
sum(sum(t,2)>0) % ncells responsive to >= 1 targ ori: 56/103

subplot(1,3,1)
imagesc(t); colorbar
title('sig>=1 is ori-tuned')
subplot(1,3,2)
imagesc(p_ttest(:,:,1)); colorbar
title('p value of small isi')
subplot(1,3,3)
imagesc(p_ttest(:,:,2)); colorbar
title('p value of large isi')
set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
saveas(gcf, ['ori_tuned_cells_across_isi.jpg'])
close

%% test no-adapter trial #

adapt_list = unique(adapterCon);
nadapt = length(adapt_list);
ntrial_cond = [];
for iadapt =  1 : nadapt % ntrial per isi is equal but not distributed evenly to every delta
    id_adapt = find(adapterCon == adapt_list(iadapt));
    
for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta)); 
    idx = intersect(id_adapt, id_delta);
    ntrial_cond(idelta, iadapt) = length(idx); % condition = specific isi (targ onset frame#) & targ dir

end
end

targ_start_list = unique(targ_start);
ngap = length(targ_start_list);
for igap =  1 : ngap % ntrial per isi is equal but not distributed evenly to every delta
    id_targ = find(targ_start == targ_start_list(igap));
    range_adapt_base = [targ_start_list(igap) - targ_stim_len : targ_start_list(igap) - 1]; % adapted baseline just bef targ onset
    range_targ_resp = [targ_start_list(igap) : targ_start_list(igap) + targ_stim_len - 1]; % targ onset til targ fin

for idelta = 1 : ndelta % ntrial per delta is equal
    id_delta = find(delta_seq == delta_list(idelta)); 
    idx = intersect(id_targ, id_delta);
    ntrial_cond2(idelta, igap) = length(idx);
end
end

ntrial_cond
ntrial_cond2

%% orientation tuning plot of indiv cell

for icell = 1 : ncell
    figure('units','normalized','outerposition',[0 0 1 1]);
    errorbar(delta_list, dfof_avg(icell,:), dfof_ste(icell,:), 'LineWidth',1)
    hold on
    if sum(sig_ttest(icell,:)) > 0
        sig_idx = sig_ttest(icell,:) > 0;
        sig_delta = delta_list(sig_idx);
        sig_star_height = dfof_avg(icell,sig_idx) + dfof_ste(icell,sig_idx) + 0.01;
        scatter(sig_delta, sig_star_height, '*', 'LineWidth',1)
    else
        sig_delta = [];
    end
    xlim([0-5, 180])
    line([0-5, 180], [0, 0], 'Color', 'g', 'LineWidth', 1);
    title(['cell ', num2str(icell), ': sensitive to ' num2str(length(sig_delta)), ' orientations'])
    saveas(gcf, ['ori_tuning_', num2str(icell)], 'jpg')
    close
end

%% ntrial_ori is similar across ori, why sig differ?

ntrial_ndelta = [];
for idelta = 1 : ndelta
    idx = find(delta_seq == delta_list(idelta)); 
    ntrial_ndelta(idelta) = length(idx); 
end

%% fit von Mises function
fit_param = pi * ones(ncell, 7);

for icell = 1 : ncell
%     figure('units','normalized','outerposition',[0 0 1 1]);
    errorbar(delta_list, dfof_avg(icell,:), dfof_ste(icell,:), 'LineWidth',1)
    hold on
    if sum(sig_ttest(icell,:)) > 0
        sig_idx = sig_ttest(icell,:) > 0;
        sig_delta = delta_list(sig_idx);
        sig_star_height = dfof_avg(icell,sig_idx) + dfof_ste(icell,sig_idx) + 0.01;
        scatter(sig_delta, sig_star_height, '*', 'LineWidth',1)
    else
        sig_delta = [];
    end
    xlim([0-5, 180+5])
    line([0-5, 180+5], [0, 0], 'Color', 'g', 'LineWidth', 1);
    
    
    theta = deg2rad(delta_list);
    data = dfof_avg(icell,:); % data = resp_avg(icell, :) - base_avg(icell, :);
    [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
    fit_param(icell,:) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
%   icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2
    
    theta_finer = deg2rad(0:1:179);
    y_fit = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(theta_finer - u1_hat))-1));
    ori_pref = rad2deg(u1_hat);
%     if ori_pref < 0
%         ori_pref = ori_pref + 180;
%     elseif ori_pref > 180
%         ori_pref = ori_pref - 180;
%     end
    ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180;
    ori_pref(ori_pref >= 180) = ori_pref(ori_pref >= 180) - 180;

    plot(rad2deg(theta_finer), y_fit, 'LineWidth', 1)
    yl = ylim; % query [ymin ymax]
    line([ori_pref, ori_pref], [yl(1), (b_hat + R1_hat)], 'Color', 'r', 'LineWidth', 1);
    
    xlabel('ori in deg')
    ylabel('dF/F')
    if isempty(sig_delta)
        title(['cell ', num2str(icell), ' has no ori-pref'])
    else
        title(['cell ', num2str(icell), ' prefer ori ', num2str(round(ori_pref, 2))])
    end
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, ['ori_tuning_fit_', num2str(icell)], 'jpg')
    close
    
end

u1_hat_cells = fit_param(:,5);
ori_pref = rad2deg(u1_hat_cells);
ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180;
ori_pref(ori_pref >= 180) = ori_pref(ori_pref >= 180) - 180;
ori_pref_cells = ori_pref;

save ori_across_cells.mat dfof_avg dfof_ste fit_param ori_pref_cells

%% bootstrap -> goodness of fit

dfof_avg_runs = pi * ones(ncell, ndelta, nrun);
dfof_ste_runs = pi * ones(ncell, ndelta, nrun);
fit_param_runs = pi * ones(ncell, 7, nrun);
ori_pref_runs = pi * ones(ncell, nrun);

theta = deg2rad(delta_list);
nrun = 1000;

for irun = 1 : nrun
    irun

    for icell = 1 : ncell
        
        for idelta = 1 : ndelta
            idx = find(delta_seq == delta_list(idelta)); 
            ntrial_cond = length(idx);
            bootstrap_draw = round(ntrial_cond * 0.7);
            idx_run = randsample(idx, bootstrap_draw, 1); % w replacement

            base_win = squeeze(tc_trials(icell, idx_run, (nOff - win_len + 1):nOff));
            base_win = mean(base_win, 2); % avg over window -> [ntrial_ori, 1]
            resp_win = squeeze(tc_trials(icell, idx_run, (trial_len - win_len + 1):trial_len));
            resp_win = mean(resp_win, 2);

%             [sig_ttest(icell, idelta), p_ttest(icell, idelta)] = ttest(base_win, resp_win, 'alpha',0.05./(ntrial_ori - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
%             base_avg(icell, idelta) = mean(base_win); % avg over sampled trials 
%             resp_avg(icell, idelta) = mean(resp_win);
%             resp_ste(icell, idelta) = std(resp_win) / sqrt(length(resp_win));

            dfof_avg_runs(icell, idelta, irun) = mean( (resp_win - base_win) ./ mean(base_win) );
            dfof_ste_runs(icell, idelta, irun) = std( (resp_win - base_win) ./ mean(base_win) ) ./ sqrt(ntrial_cond);
        end
        
        
        data = dfof_avg_runs(icell, :, irun); 
        [b_hat, k1_hat, R1_hat, u1_hat, sse, R_square] = miaovonmisesfit_ori(theta, data);
        fit_param_runs(icell, :, irun) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
    %   icell, baseline|offset, k1 sharpness, R peak response, u1 preferred orientation, sse sum of squared error, R2

%         theta_finer = deg2rad(0:1:179);
%         y_fit = b_hat + R1_hat .* exp(k1_hat.*(cos(2.*(theta_finer - u1_hat))-1));
        ori_pref = rad2deg(u1_hat);
%     if ori_pref < 0
%         ori_pref = ori_pref + 180;
%     elseif ori_pref > 180
%         ori_pref = ori_pref - 180;
%     end
        ori_pref(ori_pref < 0) = ori_pref(ori_pref < 0) + 180;
        ori_pref(ori_pref >= 180) = ori_pref(ori_pref >= 180) - 180;
        ori_pref_runs(icell, irun) = ori_pref;
        
    end
end

save ori_across_bootstrap_runs.mat dfof_avg_runs dfof_ste_runs fit_param_runs ori_pref_runs

% % sanity check
% tt = mean(dfof_avg_runs, 3);
% 
% subplot(1,2,1)
% imagesc(dfof_avg)
% colorbar
% 
% subplot(1,2,2)
% imagesc(mean(dfof_avg_runs, 3))
% % imagesc(dfof_avg_runs(:,:,1))
% colorbar

%% evaluate ori_pref across runs -> good fit

ori_closeness = sort(abs(ori_pref_runs - ori_pref_cells), 2); % sort each row aka cell
percentile_threshold = 0.90;
percentile_idx = percentile_threshold * nrun;
ori_perc = ori_closeness(:, percentile_idx);

% sum(ori_perc<22.5) / length(ori_perc)
% histogram(ori_perc,8)

good_fit_cell = ori_perc<22.5;
sum(good_fit_cell)
sig_ori_cell = sum(sig_ttest,2)>0;
sum(sig_ori_cell)
ori_cell = good_fit_cell & sig_ori_cell;
sum(ori_cell)

%% cell distribution

xname = {'segmented'; 'sig ori-tuned'; 'good fit'; 'both'};
y = [ncell, sum(sig_ori_cell), sum(good_fit_cell), sum(ori_cell)]; 
x = 1 : length(y);
bar(x, y)
set(gca,'xticklabel', xname)

for iloc = 1:numel(y)
    text(x(iloc),y(iloc), num2str(y(iloc),'%0.2f'),...
               'HorizontalAlignment','center', 'VerticalAlignment','bottom')
end

ylim([0, max(y)+10])
ylabel('number of cells')
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'number of cells.jpg')
close

%% ori tuning param dist
% fit_param(icell,:) = [icell, b_hat, k1_hat, R1_hat, u1_hat, sse, R_square];
ori_pref_qualified = rad2deg(fit_param(ori_cell, 5));
ori_pref_qualified(ori_pref_qualified<0) = ori_pref_qualified(ori_pref_qualified<0) + 180;
ori_sharp_qualified = fit_param(ori_cell, 3);

subplot(1,2,1)
edges = [delta_list, 180];
histogram(ori_pref_qualified, edges)
xlabel('preferred ori')

subplot(1,2,2)
histogram(ori_sharp_qualified, 10)
xlabel('sharpness of tuning')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'ori_pref_w_sharp.jpg')
close

%% adaptation index

