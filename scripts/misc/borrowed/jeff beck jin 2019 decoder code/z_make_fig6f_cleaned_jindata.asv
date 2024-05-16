%%

close all; clc; 
clear; clear global

root_path = 'C:\Users\ll357\Documents\inter';
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');


%% Data Preparation
% List of file names
filename = {'170323_i689_runs-002-003', '170323_i696_runs-002-003', '170324_i674_runs-002-003', ...
            '170327_i684_runs-002-003', '170503_i711_runs-002-003', '170503_i712_runs-002-003', ...
            '170510_i574_runs-002-003', '170808_i720_runs-002-003', '170810_i738_runs-002-003', ...
            '170811_i739_runs-002-003', '170816_i745_runs-002-003', '170826_i746_runs-002-003'};
        
% _newFits
for n = 1:length(filename)
    filename{n} = ['Z:\All_Staff\home\lan\Analysis\optimal decoder from jeff beck dropbox\Lindsay\', ...
        filename{n}, '_newFits.mat'];
end

PCmax = 15;
train = 3;

% Initialize data structures
DVAll.dataset = [];
DVAll.Y = [];
DVAll.cond = [];
DVAll.PV = [];

% check if session has any good units
exclude_sess = [];
ncell_good_thresh = 0; % exclude sessions with only 0-n well fit cells
theta90_wellfit_thresh = 22.5;
    
sess_kept = [];
for n = 1:length(filename)
    load([filename{n}]);
    
    [~, loc] = max(ori_fit);
    prefs{n} = loc / 180 * 2 * pi;
    idxn{n} = find(theta_90 < theta90_wellfit_thresh);
    % idxn{n} = find(well_max > 0);
    
    disp(['Dataset ', num2str(n), ' has ', num2str(length(idxn{n})), ' good units']);
    if (length(idxn{n}) <= ncell_good_thresh) % || (length(idxn{n}) > 20)
        exclude_sess = [exclude_sess, n];
    end
    sess_kept = [sess_kept, ~ismember(n, exclude_sess)];
end
sess_kept = logical(sess_kept);
filename = filename(1, sess_kept);
disp(['kept ', num2str(sum(sess_kept)), ' out of ', num2str(length(sess_kept))])


% %% construct DVAll, no rerun unless changed
% 
% % Loop through files
% for n = 1 : length(filename)
%     load([filename{n}]);
% 
%     [~, loc] = max(ori_fit);
%     prefs{n} = loc / 180 * 2 * pi;
%     idxn{n} = find(theta_90 < theta90_wellfit_thresh);
%     % idxn{n} = find(well_max > 0);
% 
%     disp(['Dataset ', num2str(n), ... 
%         ' out of ', num2str(length(filename)), ...
%         ' has ', num2str(length(idxn{n})), ' good units']);
% 
%     prefs{n} = prefs{n}(idxn{n});
% 
%     % ori_fit_empty = zeros(size(ori_fit));
%     % f{n} = ori_fit_empty(:, idxn{n}); % try empty ori_fit to check if its useless
%     f{n} = ori_fit(:, idxn{n});
% 
%     kappa{n} = abs(fft(log(f{n})));
%     kappa{n} = kappa{n}(end, :);
% 
%     % Data preparation for analysis
%     for j = 1:max(train, 2)
%         data{j}.X = [];
%         data{j}.Y = [];
% 
%         for k = 1:8
%             data{j}.X = [data{j}.X; ppResp{j, k}'];
%             data{j}.Y = [data{j}.Y; k * ones(size(ppResp{j, k}, 2), 1)];
%         end
% 
%         % Remove rows with NaN values
%         idx = ~any(isnan(data{j}.X), 2);
%         data{j}.X = data{j}.X(idx, :);
%         data{j}.Y = data{j}.Y(idx, 1);
% 
%         data{j}.X = data{j}.X(:, idxn{n});
%         data{j}.Xraw = [data{j}.X];
%     end
% 
%     % Calculate principal components
%     CC = cov([data{1}.X; data{2}.X]);
%     [V, D] = eig(CC, 'vector');
%     PCs = min(PCmax, size(data{j}.X, 2));
%     D = D(max(size(data{1}.X, 2) - PCs + 1, 1):end);
%     V = V(:, max(size(data{1}.X, 2) - PCs + 1, 1):end);
% 
%     dataAll.X = [];
%     dataAll.Xraw = [];
%     dataAll.Y = [];
%     dataAll.cond = [];
%     for j = 1:max(train, 2)
%         data{j}.X = data{j}.X * V * diag(1 ./ sqrt(D));
%         data{j}.X = data{j}.X(:, max(size(data{j}.X, 2) - PCs + 1, 1):end);
%         data{j}.X = [data{j}.X, ones(size(data{j}.X, 1), 1)];
%         dataAll.X = [dataAll.X; data{j}.X];
%         dataAll.Xraw = [dataAll.Xraw; data{j}.Xraw];
%         dataAll.Y = [dataAll.Y; data{j}.Y];
%         dataAll.cond = [dataAll.cond; j * ones(size(data{j}.Y))];
%     end
% 
%     % Calculate DVs
%     % [~, DVAlltemp, Btemp, ~] = getDVs(data, dataAll, 0, prefs{n}(:,:), kappa{n}(:,:), f{n}(:,:));
%     [DVtemp, DVAlltemp, Btemp, ~] = getDVs(data, dataAll, 0, prefs{n}(:,:), kappa{n}(:,:), f{n}(:,:)); % like alldataanalysis.m
% 
%     for j = 1:max(train, 2)
%         Btemp{j} = Btemp{j}(1:end - 1, 1);
%     end
% 
%     % Remove bias from estimators
%     for j = 1:max(train, 2)
%         DVAll.dataset = [DVAll.dataset; n * ones(size(DVAlltemp{j}.opt))];
%         DVAll.Y = [DVAll.Y; data{j}.Y];
%         DVAll.cond = [DVAll.cond; j * ones(size(data{j}.Y));];
%         DVAll.PV = [DVAll.PV; DVAlltemp{j}.PV];
%         % DVAll.PV = [DVAll.PV; DVtemp{j}.PVemp];
%     end
% end

%% load dvall

cd('C:\Users\ll357\Documents\inter\results\decoder_grat8\pop vec decoder jin2019 jeff')
tmp = load('pop_vec_decoder_jin2019_dvall.mat');
DVAll = tmp.DVAll;

%% ref0 decoder (default: across ISI)
% % decoder ori=0 vs other

ori_fp = [1, 2, 3, 4, 5, 6, 7]; % TODO: why??
NDC = 500;
dv = [0:NDC] / NDC;

clear AUROC
kk = 0;
k = 8; % leftover variable from DVAll construction, not sure why it is necessary
usedatasets=[2:4,6:12]; % original dataset filter by jeff
for dataset = usedatasets
    kk = kk + 1;
    for j = 1:3
        for difficulty = 1:5
            idxfp = (DVAll.Y == 8 & DVAll.cond == 1 & logical(sum(DVAll.dataset == dataset, 2)));
            switch difficulty
                case 1
                    not8 = 8;
                case 2
                    not8 = [1, 7];
                case 3
                    not8 = [2, 6];
                case 4
                    not8 = [3, 5];
                case 5
                    not8 = 4;
            end
            idxcd = (logical(sum(DVAll.Y == not8, 2)) & DVAll.cond == j & ...
                logical(sum(DVAll.dataset == dataset, 2)));
            DVtemp = abs(getfield(DVAll, 'PV'));
            CD = mean(DVtemp(idxcd) > dv);
            FP = mean(DVtemp(idxfp) >= dv);

            AUROC{k}(kk, j, difficulty) = -trapz(FP, CD);
        end
    end
end
AUROC_ref0_acrossISI = AUROC;

%% ref0 decoder (modified: within ISI)

ori_fp = [1, 2, 3, 4, 5, 6, 7]; % TODO: why??
NDC = 500;
dv = [0:NDC] / NDC;

clear AUROC
kk = 0;
k = 8; % leftover variable from DVAll construction, not sure why it is necessary
usedatasets=[2:4,6:12]; % original dataset filter by jeff
for dataset = usedatasets
    kk = kk + 1;
    for j = 1:3
        for difficulty = 1:5
            idxfp = (DVAll.Y == 8 & DVAll.cond == j & logical(sum(DVAll.dataset == dataset, 2)));
            switch difficulty
                case 1
                    not8 = 8;
                case 2
                    not8 = [1, 7];
                case 3
                    not8 = [2, 6];
                case 4
                    not8 = [3, 5];
                case 5
                    not8 = 4;
            end
            idxcd = (logical(sum(DVAll.Y == not8, 2)) & DVAll.cond == j & ...
                logical(sum(DVAll.dataset == dataset, 2)));
            DVtemp = abs(getfield(DVAll, 'PV'));
            CD = mean(DVtemp(idxcd) > dv);
            FP = mean(DVtemp(idxfp) >= dv);

            AUROC{k}(kk, j, difficulty) = -trapz(FP, CD);
        end
    end
end
AUROC_ref0_withinISI = AUROC;

%% neighbor decoder (default: within ISI)
% % decode each ori against its left neighbor (sorted by ori_dist from adapter)
% % when ori=0, decode against itself

NDC = 500;
dv = [0:NDC] / NDC;

clear AUROC
kk=0;
usedatasets=[2:4,6:12]; % original dataset filter by jeff
for dataset = usedatasets
    kk = kk + 1;
    for j = 1:3
        for pair = 1:7 % despite ori_dist having 5 cases, we need to separate neighboring pairs (2-1 vs 6-7)
            switch pair
                case 1
                    ori_cd = 8;
                    ori_fp = 8;
                case 2
                    ori_cd = [1, 7];
                    ori_fp = 8;
                case 3
                    ori_cd = 2; % NOTE: fold up case 3 & 4
                    ori_fp = 1;
                case 4
                    ori_cd = 6;
                    ori_fp = 7;
                case 5
                    ori_cd = 3; % NOTE: fold up case 5 & 6
                    ori_fp = 2;
                case 6
                    ori_cd = 5;
                    ori_fp = 6;
                case 7
                    ori_cd = 4;
                    ori_fp = [3, 5];
            end
            idxfp = (logical(sum(DVAll.Y == ori_fp, 2)) ...
                & DVAll.cond == j ... % NOTE: in neighbor task, compare within isi condition, instead of always comparing to isi=250 ori=0
                & logical(sum(DVAll.dataset == dataset, 2)));

            idxcd = (logical(sum(DVAll.Y == ori_cd, 2)) ...
                & DVAll.cond == j ...
                & logical(sum(DVAll.dataset == dataset, 2)));

            DVtemp = abs(getfield(DVAll, 'PV'));
            CD = mean(DVtemp(idxcd) > dv);
            FP = mean(DVtemp(idxfp) >= dv);

            AUROC{k}(kk, j, pair) = -trapz(FP, CD);
        end
    end
end
AUROC_neighbor_withinISI = AUROC;

%% neighbor decoder (modified: across ISI)

NDC = 500;
dv = [0:NDC] / NDC;

clear AUROC
kk=0;
usedatasets=[2:4,6:12]; % original dataset filter by jeff
for dataset = usedatasets
    kk = kk + 1;
    for j = 1:3
        for pair = 1:7 % despite ori_dist having 5 cases, we need to separate neighboring pairs (2-1 vs 6-7)
            switch pair
                case 1
                    ori_cd = 8;
                    ori_fp = 8;
                case 2
                    ori_cd = [1, 7];
                    ori_fp = 8;
                case 3
                    ori_cd = 2; % NOTE: fold up case 3 & 4
                    ori_fp = 1;
                case 4
                    ori_cd = 6;
                    ori_fp = 7;
                case 5
                    ori_cd = 3; % NOTE: fold up case 5 & 6
                    ori_fp = 2;
                case 6
                    ori_cd = 5;
                    ori_fp = 6;
                case 7
                    ori_cd = 4;
                    ori_fp = [3, 5];
            end
            idxfp = (logical(sum(DVAll.Y == ori_fp, 2)) ...
                & DVAll.cond == 1 ... % NOTE: across ISI
                & logical(sum(DVAll.dataset == dataset, 2)));

            idxcd = (logical(sum(DVAll.Y == ori_cd, 2)) ...
                & DVAll.cond == j ...
                & logical(sum(DVAll.dataset == dataset, 2)));

            DVtemp = abs(getfield(DVAll, 'PV'));
            CD = mean(DVtemp(idxcd) > dv);
            FP = mean(DVtemp(idxfp) >= dv);

            AUROC{k}(kk, j, pair) = -trapz(FP, CD);
        end
    end
end
AUROC_neighbor_acrossISI = AUROC;

%% stats & plot: within ISI

close all
for decoder_mode = 0:1
    if decoder_mode == 0
        AUROC = AUROC_ref0_withinISI;
        decoder_mode_str = 'ref0 withinISI'
    elseif decoder_mode == 1
        AUROC = AUROC_neighbor_withinISI;
        decoder_mode_str = 'neighbor withinISI'
    end
    
    tmp_250 = squeeze(AUROC{k}(:, 1, :));
    tmp_750 = squeeze(AUROC{k}(:, 2, :));
    tmp_6000 = squeeze(AUROC{k}(:, 3, :)); % ISI order for jin data: 250-750-6000

    % tmp_250 = squeeze(AUROC{k}(:, 1, :));
    % tmp_6000 = squeeze(AUROC{k}(:, 2, :)); % ISI order for my data: 250-6000-750
    % tmp_750 = squeeze(AUROC{k}(:, 3, :)); % determined by dfof_resp_trialwise_jeff.m
    
    if decoder_mode == 1
        tmp_250_fold = [tmp_250(:, 1:2), mean(tmp_250(:, 3:4), 2), mean(tmp_250(:, 5:6), 2), tmp_250(:, end)];
        tmp_750_fold = [tmp_750(:, 1:2), mean(tmp_750(:, 3:4), 2), mean(tmp_750(:, 5:6), 2), tmp_750(:, end)];
        tmp_6000_fold = [tmp_6000(:, 1:2), mean(tmp_6000(:, 3:4), 2), mean(tmp_6000(:, 5:6), 2), tmp_6000(:, end)];

        tmp_250 = tmp_250_fold;
        tmp_750 = tmp_750_fold;
        tmp_6000 = tmp_6000_fold;
    end
    
    [~, p750, ~, ~] = ttest(tmp_250(:, 2), tmp_750(:, 2), "Tail","right");
    disp([decoder_mode_str, ' decoder 22-0 250 vs 750 sig ', num2str(p750)])
    [~, p6000, ~, ~] = ttest(tmp_250(:, 2), tmp_6000(:, 2), "Tail","right");
    disp([decoder_mode_str, ' decoder 22-0 250 vs 6000 sig ', num2str(p6000)])
    

    xa = [0, 22.5, 45, 67.5, 90];
    figure
    subplot(121)
    errorbar(xa, nanmean(tmp_250), ...
                nanstd(tmp_250) / sqrt(size(tmp_250, 1)), 'b')
    hold on
    errorbar(xa+1, nanmean(tmp_750), ...
                nanstd(tmp_750) / sqrt(size(tmp_250, 1)), 'r')
    ylabel('AUROC')
    xlabel('Orientation difference')
    axis([-5, 95, 0.4, 1])
    legend('250', '750', 'Location','northwest')
    title([decoder_mode_str, ' decoder'])


    subplot(122)
    errorbar(xa, nanmean(tmp_250), ...
                nanstd(tmp_250) / sqrt(size(tmp_250, 1)), 'b')
    hold on
    errorbar(xa+1, nanmean(tmp_6000), ...
                nanstd(tmp_6000) / sqrt(size(tmp_6000, 1)), 'r')
    ylabel('AUROC')
    xlabel('Orientation difference')
    axis([-5, 95, 0.4, 1])
    legend('250', '6000', 'Location','northwest')
    title([decoder_mode_str, ' decoder'])
  
end


%% stats & plot: across ISI (false positive trial index = cond isi 250)

for decoder_mode = 0:1
    if decoder_mode == 0
        AUROC = AUROC_ref0_acrossISI;
        decoder_mode_str = 'ref0 acrossISI'
    elseif decoder_mode == 1
        AUROC = AUROC_neighbor_acrossISI;
        decoder_mode_str = 'neighbor acrossISI'
    end
    
    tmp_250 = squeeze(AUROC{k}(:, 1, :));
    tmp_750 = squeeze(AUROC{k}(:, 2, :));
    tmp_6000 = squeeze(AUROC{k}(:, 3, :)); % ISI order for jin data: 250-750-6000

    % tmp_250 = squeeze(AUROC{k}(:, 1, :));
    % tmp_6000 = squeeze(AUROC{k}(:, 2, :)); % ISI order for my data: 250-6000-750
    % tmp_750 = squeeze(AUROC{k}(:, 3, :)); % determined by dfof_resp_trialwise_jeff.m
    
    if decoder_mode == 1
        tmp_250_fold = [tmp_250(:, 1:2), mean(tmp_250(:, 3:4), 2), mean(tmp_250(:, 5:6), 2), tmp_250(:, end)];
        tmp_750_fold = [tmp_750(:, 1:2), mean(tmp_750(:, 3:4), 2), mean(tmp_750(:, 5:6), 2), tmp_750(:, end)];
        tmp_6000_fold = [tmp_6000(:, 1:2), mean(tmp_6000(:, 3:4), 2), mean(tmp_6000(:, 5:6), 2), tmp_6000(:, end)];

        tmp_250 = tmp_250_fold;
        tmp_750 = tmp_750_fold;
        tmp_6000 = tmp_6000_fold;
    end
    
    [~, p750, ~, ~] = ttest(tmp_250(:, 2), tmp_750(:, 2), "Tail","right");
    disp([decoder_mode_str, ' decoder 22-0 250 vs 750 sig ', num2str(p750)])
    [~, p6000, ~, ~] = ttest(tmp_250(:, 2), tmp_6000(:, 2), "Tail","right");
    disp([decoder_mode_str, ' decoder 22-0 250 vs 6000 sig ', num2str(p6000)])
    

    xa = [0, 22.5, 45, 67.5, 90];
    figure
    subplot(121)
    errorbar(xa, nanmean(tmp_250), ...
                nanstd(tmp_250) / sqrt(size(tmp_250, 1)), 'b')
    hold on
    errorbar(xa+1, nanmean(tmp_750), ...
                nanstd(tmp_750) / sqrt(size(tmp_250, 1)), 'r')
    ylabel('AUROC')
    xlabel('Orientation difference')
    axis([-5, 95, 0.4, 1])
    legend('250', '750', 'Location','northwest')
    title([decoder_mode_str, ' decoder'])


    subplot(122)
    errorbar(xa, nanmean(tmp_250), ...
                nanstd(tmp_250) / sqrt(size(tmp_250, 1)), 'b')
    hold on
    errorbar(xa+1, nanmean(tmp_6000), ...
                nanstd(tmp_6000) / sqrt(size(tmp_6000, 1)), 'r')
    ylabel('AUROC')
    xlabel('Orientation difference')
    axis([-5, 95, 0.4, 1])
    legend('250', '6000', 'Location','northwest')
    title([decoder_mode_str, ' decoder'])
  
end
