%% same img session as in tutorial

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

%% load data (took the ref TC)

fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName); % load behavior data, aka "input"

CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata, aka "info" % check content

tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs.mat']);

%% use npSub_tc for further analysis

nOn = input.nScansOn; 
nOff = input.nScansOff;
trial_len = nOn + nOff;
ntrials = size(input.tGratingDirectionDeg,2);
size(npSub_tc) % nframe * ncell
ncell = size(npSub_tc, 2);

tc_trials = reshape(npSub_tc,[nOn+nOff, ntrials, ncell]); 
tc_trials = permute(tc_trials, [3,2,1]);
size(tc_trials) % ncell * ntrial * trial_len

Dir = celleqel2mat_padded(input.tGratingDirectionDeg); 
convert_idx = Dir>=180;
Ori = Dir;
Ori(convert_idx) = Ori(convert_idx) - 180; 
Ori_list = unique(Ori);
nOri = length(Ori_list);

%% probe win_len

tt = mean(mean(tc_trials,1),2);
plot(1:trial_len, squeeze(tt)) 
% trial = 0-60 off + 61-90 on. signal decays >2/3 after 0-30 off
% as short as 3 frames would suffice. now take 10 frames as window len: 51-60 = base, 81-90 = resp
win_len = 10;

%% cells sensitive to orientations
sig_ttest = pi * ones(ncell, nOri);
p_ttest = pi * ones(ncell, nOri);
base_avg = pi * ones(ncell, nOri);
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
        
        [sig_ttest(icell, iOri), p_ttest(icell, iOri)] = ttest(base_win, resp_win, 'alpha',0.05./(ntrials_ori - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
        base_avg(icell, iOri) = mean(base_win); % avg over trials of same ori
        resp_avg(icell, iOri) = mean(resp_win);
        resp_ste(icell, iOri) = std(resp_win) / sqrt(length(resp_win));

        dfof_avg(icell, iOri) = mean( (resp_win - base_win) ./ mean(base_win) );
        dfof_ste(icell, iOri) = std( (resp_win - base_win) ./ mean(base_win) ) / sqrt(ntrials_ori);
    end
end

sum(sum(sig_ttest,2)>0) % ncells responsive to >= 1 orientation: 80->89/148 after set tail

base = mean(tc_trials(:,:, (nOff - win_len + 1):nOff), 3);
resp = mean(tc_trials(:,:, (trial_len - win_len + 1):trial_len), 3);
df = resp - base;

%% orientation tuning plot of indiv cell

for icell = 1 : ncell
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
    line([0-5, 180], [0, 0], 'Color', 'g', 'LineWidth', 1);
    title(['cell ', num2str(icell), ': sensitive to ' num2str(length(sig_ori)), ' orientations'])
    saveas(gcf, ['ori_tuning_', num2str(icell)], 'jpg')
    close
end

%% ntrial_ori is similar across ori, why sig differ?

ntrials_nOri = [];
for iOri = 1 : nOri
    idx = find(Ori == Ori_list(iOri)); 
    ntrials_nOri(iOri) = length(idx); 
end

%% fit von Mises function
fit_param = pi * ones(ncell, 7);

for icell = 1 : ncell
%     figure('units','normalized','outerposition',[0 0 1 1]);
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
    xlim([0-5, 180+5])
    line([0-5, 180+5], [0, 0], 'Color', 'g', 'LineWidth', 1);
    
    
    theta = deg2rad(Ori_list);
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
    if isempty(sig_ori)
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

dfof_avg_runs = pi * ones(ncell, nOri, nrun);
dfof_ste_runs = pi * ones(ncell, nOri, nrun);
fit_param_runs = pi * ones(ncell, 7, nrun);
ori_pref_runs = pi * ones(ncell, nrun);

theta = deg2rad(Ori_list);
nrun = 1000;

for irun = 1 : nrun
    irun

    for icell = 1 : ncell
        
        for iOri = 1 : nOri
            idx = find(Ori == Ori_list(iOri)); 
            ntrials_ori = length(idx);
            bootstrap_draw = round(ntrials_ori * 0.7);
            idx_run = randsample(idx, bootstrap_draw, 1); % w replacement

            base_win = squeeze(tc_trials(icell, idx_run, (nOff - win_len + 1):nOff));
            base_win = mean(base_win, 2); % avg over window -> [ntrial_ori, 1]
            resp_win = squeeze(tc_trials(icell, idx_run, (trial_len - win_len + 1):trial_len));
            resp_win = mean(resp_win, 2);

%             [sig_ttest(icell, iOri), p_ttest(icell, iOri)] = ttest(base_win, resp_win, 'alpha',0.05./(ntrials_ori - 1), 'tail', 'left'); % sig = base<resp, Bonferroni correction
%             base_avg(icell, iOri) = mean(base_win); % avg over sampled trials 
%             resp_avg(icell, iOri) = mean(resp_win);
%             resp_ste(icell, iOri) = std(resp_win) / sqrt(length(resp_win));

            dfof_avg_runs(icell, iOri, irun) = mean( (resp_win - base_win) ./ mean(base_win) );
            dfof_ste_runs(icell, iOri, irun) = std( (resp_win - base_win) ./ mean(base_win) ) ./ sqrt(ntrials_ori);
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

save ori_bootstrap.mat dfof_avg_runs dfof_ste_runs fit_param_runs ori_pref_runs

% san
tt = mean(dfof_avg_runs, 3);

subplot(1,2,1)
imagesc(dfof_avg)
colorbar

subplot(1,2,2)
imagesc(mean(dfof_avg_runs, 3))
% imagesc(dfof_avg_runs(:,:,1))
colorbar

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
edges = [Ori_list, 180];
histogram(ori_pref_qualified, edges)
xlabel('preferred ori')

subplot(1,2,2)
histogram(ori_sharp_qualified, 10)
xlabel('sharpness of tuning')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, 'ori_pref_w_sharp.jpg')
close
