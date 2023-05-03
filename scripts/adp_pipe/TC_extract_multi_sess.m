%% readme
% this script register a multi-session dataset, 
% then extract timecourse from cells segmented by cellpose
% it only support ONE single day multi sess at the moment
% refactored to loop over days: 2023-05-03

clear all
clearvars –global
close all
clc

% write text output from command window to log file
cd('C:\Users\ll357\Documents\inter\data\')
diary matlab_log % batch cellpose mix14 (grat_SF6_allen_nat8)
datetime

%% setup

dir_meta = 'Z:\All_Staff\home\lan\Data\2P_images\mat_inter/adp_dataset_master.xlsx';
dataset_meta = readtable(dir_meta);

% clearvars -except dataset_meta
% clearvars –global

stim_type = 'bunny' % grat_SF6_allen_nat8 will pretend to be bunny, since they share bunnies6.mwel
dataset_table = dataset_meta(strcmp(dataset_meta.paradigm, 'grat_SF6_allen_nat8'), :);
% dataset_table = dataset_meta(dataset_meta.date == 230330, :)

ndate = length(unique(dataset_table.date));
date_arr = unique(dataset_table.date);

for idate = 1:ndate

dataset_date  = dataset_table(dataset_table.date == date_arr(idate), :);
nsess = size(dataset_date, 1);

%% draw cellpose_tif for each sess in each date

for isess = 1:nsess

isess, nsess
dataset_now = dataset_date(isess,:);
arg_mouse = dataset_now.mouse
arg_date = num2str(dataset_now.date)
arg_ImgFolder = dataset_now.num{1};

ref_ImgFolder = dataset_date(1,:).num{1}
run_str_ref = ['runs-' ref_ImgFolder(1,:)]
if isess == 1
    disp('take first session 002 of multisess data as ref session')
end

imouse = ['i', num2str(arg_mouse)];
dir_analysis = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse, ...
    '\', arg_date, '_', imouse, '_runs-', arg_ImgFolder];
try cd(dir_analysis)
catch
    mkdir(dir_analysis)
    cd(dir_analysis)
end

% % check if cellpose tif exists in date folder, bc TC is in sess folder
if ~isempty(dir('cellpose_stim_resp_gauss.tif'))
    disp('cellpose tif exists, skip to next set:')
    continue
end

disp('add your stim_type into function: bunny, mix, grat6, grating(grat1)')
% % register every sess against data_avg of session 002
[~, ~, date, imouse, run_str] = get_data_reg_cellpose_tif(...
    arg_mouse, arg_date, arg_ImgFolder, stim_type, run_str_ref); 
disp(['got data_reg & cellpose tif for ', ...
    num2str(arg_mouse), ' ', arg_date, ' ', arg_ImgFolder])
close all

end

%% merge into multisess tif

dir_final_tif = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse];
cd(dir_final_tif) % one level up from sess folder
file_list = dir(fullfile(dir_final_tif, '**\data_dfof.mat'));
assert(size(file_list, 1) == nsess); % ensure each session saved a mat for cellpose tiff

tmp = pi * ones(264, 796, nsess);
for i = 1 : length(file_list)
    file_name = [file_list(i).folder, '\', file_list(i).name];
    load(file_name, 'data_dfof_max')
    tmp(:,:,i) = data_dfof_max;
end

data_dfof_multisess = mean(tmp, 3); % aggregate sess tif to get final tif
save_mat_as_tif(data_dfof_multisess) % pass to cellpose in ipynb, who reads from multisess tif folder (one level above sess folder)

%% wait for cellpose, then write timecourse

disp('wait for batch_cellpose.ipynb')
while ~exist('cellpose_mask.mat','file')
    pause(60) % wait for cellpose to run in ipynb
end
disp('got cellpose mask, now extract TC from each sess')

file_list = dir(fullfile(dir_final_tif, '**\*_reg_shifts.mat'));
file_list.name
for i = 1 : length(file_list)
    if ~isempty(dir([file_list(i).folder, '\', '*TCs_cellpose.mat'])) % pass if multisess cellpose time course exist
        disp('sess cellpose time course exists, skip to next set:')
        disp(i+1)
        continue
    end

    file_name = [file_list(i).folder, '\', file_list(i).name];
    load(file_name, 'out')
    arg_ImgFolder = dataset_date(i,:).num{1};

    [data, ~, ~, ~, ~, run_str_sess] = load_sbx_data(arg_mouse, arg_date, arg_ImgFolder);
    [outs, data_reg] = stackRegister_MA_LL(double(data), [], [], out); % re-register to get data_reg back
    clear data
    
    cd(dir_final_tif)
    tif_name = [dir_final_tif, '\cellpose_mask.mat'];
    LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan'; % TODO: spare generating LL_base from function load_sbx_data
    npSub_tc = get_cellpose_timecourse(data_reg, tif_name, ...
        LL_base, arg_date, imouse, run_str_sess);
    
    disp([num2str(isess), ' sess done out of ', num2str(nsess)])
    clear data_reg
    clear global % suspect weird trace is bc residual global var affecting sbxread

end

end
