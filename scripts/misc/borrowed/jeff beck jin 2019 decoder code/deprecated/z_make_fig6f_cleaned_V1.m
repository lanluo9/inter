%%

close all; clc; 
clear; clear global

root_path = 'C:\Users\ll357\Documents\inter';
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

dir_meta = 'Z:\All_Staff\home\lan\Data\2P_images\mat_inter/adp_dataset_master.xlsx';
dataset_meta = readtable(dir_meta);

stim_type = 'grating' % grat_8ori_3isi
dataset_table = dataset_meta(strcmp(dataset_meta.paradigm, stim_type), :);

dataset_table = dataset_table(strcmp(dataset_table.gcamp, '6s'), :);
seg_bool = dataset_table.manual_seg | dataset_table.cellpose_seg; % exclude not-segmented data
dataset_table = dataset_table(seg_bool, :);

select_area = 'LI';
area_bool = logical(strcmp(dataset_table.area, select_area));
dataset_table = dataset_table(area_bool, :);

sum(strcmp(dataset_table.area, 'V1'))
sum(strcmp(dataset_table.area, 'LM'))
sum(strcmp(dataset_table.area, 'LI'))

nset = size(dataset_table, 1);
filename = cell(1, nset);
for iset = 1:nset
    dataset_now = dataset_table(iset,:);
    arg_mouse = dataset_now.mouse;
    arg_date = num2str(dataset_now.date);
    arg_ImgFolder = dataset_now.num{1};
    area = dataset_now.area{1};
    
    imouse = ['i', num2str(arg_mouse)];
    area_mouse_date_sess = [area '_' imouse '_' arg_date '_' arg_ImgFolder];
    mapped_path = 'Z:\All_Staff\home\lan\Data\2P_images';
    try
        segment_suffix = '_cellpose';
        result_folder = [mapped_path, '\mat_inter\', area_mouse_date_sess, segment_suffix];
        cd(result_folder)
    catch
        segment_suffix = '';
        result_folder = [mapped_path, '\mat_inter\', area_mouse_date_sess, segment_suffix];
        cd(result_folder)
    end
    jeff_file = fullfile(result_folder, 'pop_vec_decoder_jeff_visp.mat');
    filename{1, iset} = jeff_file;
end

%% Data Preparation
% List of file names
% filename = {'170323_i689_runs-002-003', '170323_i696_runs-002-003', '170324_i674_runs-002-003', ...
%             '170327_i684_runs-002-003', '170503_i711_runs-002-003', '170503_i712_runs-002-003', ...
%             '170510_i574_runs-002-003', '170808_i720_runs-002-003', '170810_i738_runs-002-003', ...
%             '170811_i739_runs-002-003', '170816_i745_runs-002-003', '170826_i746_runs-002-003'};
        
PCmax = 15;
train = 3;

% Initialize data structures
DVAll.dataset = [];
DVAll.Y = [];
DVAll.cond = [];
DVAll.PV = [];

% check if session has any good units
exclude_sess = [];
if strcmp(select_area, 'LM')
    exclude_sess = [18]; % too many well fit cells, unlikely to be LM
end

ncell_wellfit_thresh = 5; % exclude sessions with only 0-n well fit cells
theta90_wellfit_thresh = 22.5;
if strcmp(select_area, 'LI')
    theta90_wellfit_thresh = 45; % relax well fit criteria for LI
end
    
sess_kept = [];
for n = 1:length(filename)
    % Load data from files
    load([filename{n}]);
    
    [~, loc] = max(ori_fit);
    prefs{n} = loc / 180 * 2 * pi;
    idxn{n} = find(theta_90 < theta90_wellfit_thresh);
    
    % Display the number of good units using theta_90
    disp(['Dataset ', num2str(n), ' has ', num2str(length(idxn{n})), ' good units using theta_90']);
    if length(idxn{n}) <= ncell_wellfit_thresh
        % disp(['skipping dataset ', num2str(n), ' ', filename{n}, ' due to 0 good (well fit) units'])
        exclude_sess = [exclude_sess, n];
    end
    sess_kept = [sess_kept, ~ismember(n, exclude_sess)];
end
sess_kept = logical(sess_kept);
filename = filename(1, sess_kept);
disp(['kept ', num2str(sum(sess_kept)), ' out of ', num2str(length(sess_kept)), ...
    ' in ', select_area])

%%

% Loop through files
for n = 1:length(filename)
    % Load data from files
    load([filename{n}]);
    
    [~, loc] = max(ori_fit);
    prefs{n} = loc / 180 * 2 * pi;
    idxn{n} = find(theta_90 < theta90_wellfit_thresh);
    
    % Display the number of good units using theta_90
    disp(['Dataset ', num2str(n), ... 
        ' out of ', num2str(length(filename)), ...
        ' has ', num2str(length(idxn{n})), ' good units using theta_90']);
    
    prefs{n} = prefs{n}(idxn{n});
    f{n} = ori_fit(:, idxn{n});
    kappa{n} = abs(fft(log(f{n})));
    kappa{n} = kappa{n}(end, :);

    % Data preparation for analysis
    for j = 1:max(train, 2)
        data{j}.X = [];
        data{j}.Y = [];

        for k = 1:8
            data{j}.X = [data{j}.X; ppResp{j, k}'];
            data{j}.Y = [data{j}.Y; k * ones(size(ppResp{j, k}, 2), 1)];
        end

        % Remove rows with NaN values
        idx = ~any(isnan(data{j}.X), 2);
        data{j}.X = data{j}.X(idx, :);
        data{j}.Y = data{j}.Y(idx, 1);

        data{j}.X = data{j}.X(:, idxn{n});
        data{j}.Xraw = [data{j}.X];
    end

    % Calculate principal components
    CC = cov([data{1}.X; data{2}.X]);
    [V, D] = eig(CC, 'vector');
    PCs = min(PCmax, size(data{j}.X, 2));
    D = D(max(size(data{1}.X, 2) - PCs + 1, 1):end);
    V = V(:, max(size(data{1}.X, 2) - PCs + 1, 1):end);

    dataAll.X = [];
    dataAll.Xraw = [];
    dataAll.Y = [];
    dataAll.cond = [];
    for j = 1:max(train, 2)
        data{j}.X = data{j}.X * V * diag(1 ./ sqrt(D));
        data{j}.X = data{j}.X(:, max(size(data{j}.X, 2) - PCs + 1, 1):end);
        data{j}.X = [data{j}.X, ones(size(data{j}.X, 1), 1)];
        dataAll.X = [dataAll.X; data{j}.X];
        dataAll.Xraw = [dataAll.Xraw; data{j}.Xraw];
        dataAll.Y = [dataAll.Y; data{j}.Y];
        dataAll.cond = [dataAll.cond; j * ones(size(data{j}.Y))];
    end

    % Calculate DVs
    [~, DVAlltemp, Btemp, ~] = getDVs(data, dataAll, 0, prefs{n}(:,:), kappa{n}(:,:), f{n}(:,:));

    for j = 1:max(train, 2)
        Btemp{j} = Btemp{j}(1:end - 1, 1);
    end

    % Remove bias from estimators
    for j = 1:max(train, 2)
        DVAll.dataset = [DVAll.dataset; n * ones(size(DVAlltemp{j}.opt))];
        DVAll.Y = [DVAll.Y; data{j}.Y];
        DVAll.cond = [DVAll.cond; j * ones(size(data{j}.Y));];
        DVAll.PV = [DVAll.PV; DVAlltemp{j}.PV];
    end
end

%% auroc
not8 = [1, 2, 3, 4, 5, 6, 7]; % TODO: why??
NDC = 500;
dv = [0:NDC] / NDC;
usedatasets = [1 : max(DVAll.dataset)];
norm_ndata = max(DVAll.dataset);
kk = 0;

clear AUROC
% for k = 1:8
for dataset = usedatasets
    kk = kk + 1;
    for j = 1:2
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
% end

%% Plotting

% tmp = load('pop_vec_decoder_jeff_res_LI_visp.mat');
% AUROC = tmp.AUROC;
% norm_ndata = size(tmp.AUROC{1, 8}, 1);
% % norm_ndata = sqrt(norm_ndata);

xa = [0, 22.5, 45, 67.5, 90];
tmp_250 = squeeze(AUROC{k}(:, 1, :));
tmp_750 = squeeze(AUROC{k}(:, 2, :));

% close all
figure;
hold on

errorbar(xa, nanmean(tmp_250), ...
            nanstd(tmp_250) / norm_ndata, 'b')
errorbar(xa+1, nanmean(tmp_750), ...
            nanstd(tmp_750) / norm_ndata, 'r')

% errorbar(xa, nanmedian(tmp_250), ...
%             nanstd(tmp_250) / norm_ndata, 'b')
% errorbar(xa+1, nanmedian(tmp_750), ...
%             nanstd(tmp_750) / norm_ndata, 'r')

title('PV')
ylabel('AUROC')
xlabel('Orientation difference')
axis([-5, 95, 0.4, 1])
legend('250', '750', 'Location','southeast')

cd('C:\Users\ll357\Documents\inter\results\decoder_grat8\pop vec decoder jin2019 jeff')
% save pop_vec_decoder_jeff_res_LM_visp.mat AUROC

%%
