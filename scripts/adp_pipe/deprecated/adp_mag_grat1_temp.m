
clear all
close all
clc

%% setup

% database_path = 'C:/Users/GlickfeldLab/Documents/test/inter/';
database_path = 'C:\Users\ll357\Documents\inter\';
master_xls = [database_path, 'data/mix50_grat1.csv'];
dataset_meta = readtable(master_xls);

clearvars -except dataset_meta
clearvars –global

stim_type = 'grat6' % only 1 grating, but modified from twoStim 2p frames mwel, just like grat SF6

% dataset_mix = dataset_meta(ismember(dataset_meta.stim_type, 'mix'),:);
% dataset_grat1 = dataset_meta(ismember(dataset_meta.stim_type, 'grating'),:);
% dataset_grat6 = dataset_meta(ismember(dataset_meta.stim_type, 'grat6'),:);

% if strcmp(stim_type, 'mix')
%     dataset_table = dataset_mix;
% elseif strcmp(stim_type, 'grating')
%     dataset_table = dataset_grat1
% elseif strcmp(stim_type, 'grat6')
%     dataset_table = dataset_grat6
% end

dataset_table = dataset_meta;
dataset_table = dataset_table(dataset_table.date == 230207, :)
nset = size(dataset_table,1);

%% draw cellpose_tif for each sess

% TODO: refactor to read mwel name from "archive 2P Imaging Notes Lan" -> 
% modify get_data_reg_cellpose_tif.m to auto fill Unrecognized field name "cStimOn" etc

for iset = 1:nset

iset, nset
dataset_now = dataset_table(iset,:);
arg_mouse = dataset_now.mouse
arg_date = num2str(dataset_now.date)
arg_ImgFolder = dataset_now.num(1);
arg_ImgFolder = ['00', num2str(arg_ImgFolder)]
ref_ImgFolder = ['00', num2str(dataset_table(1,:).num(1))];
run_str_ref = ['runs-' ref_ImgFolder(1,:)]

if iset == 1 % strcmp(stim_type, 'mix') && 
    disp('take first session 002 of mix50 as ref session')
end

imouse = ['i', num2str(arg_mouse)];
dir_analysis = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse, ...
    '\', arg_date, '_', imouse, '_runs-', arg_ImgFolder];
try cd(dir_analysis)
catch
    mkdir(dir_analysis)
    cd(dir_analysis)
end

if ~isempty(dir('cellpose_stim_resp_gauss.tif'))
    disp('cellpose tif exists, skip to next set:')
    continue
end

disp('add your stim_type into function: bunny, mix, grat6, grating(grat1)')
[~, ~, date, imouse, run_str] = get_data_reg_cellpose_tif(... % register every sess against data_avg of session 002
    arg_mouse, arg_date, arg_ImgFolder, stim_type, run_str_ref); 
disp(['got data_reg & cellpose tif for session ', arg_ImgFolder])

end

%% merge into multisess tif

dir_final_tif = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse];
cd(dir_final_tif) % one level up from sess folder
file_list = dir(fullfile(dir_final_tif, '**\data_dfof.mat'));
assert(size(file_list, 1) == nset); % ensure each session saved a mat for cellpose tiff

tmp = pi * ones(264, 796, nset);
for i = 1 : length(file_list)
    file_name = [file_list(i).folder, '\', file_list(i).name];
    load(file_name, 'data_dfof_max')
    tmp(:,:,i) = data_dfof_max;
end

data_dfof_multisess = mean(tmp, 3); % aggregate sess tif to get final tif
save_mat_as_tif(data_dfof_multisess) % pass to cellpose in ipynb, who reads from multisess tif folder (one level above sess folder)

%%

while ~exist('cellpose_mask.mat','file')
    pause(60) % wait for cellpose to run in ipynb
    toc
end
disp('got cellpose mask, now extract TC from each sess')

file_list = dir(fullfile(dir_final_tif, '**\*_reg_shifts.mat'));
file_list.name
for i = 1 : (length(file_list))
    if ~isempty(dir([file_list(i).folder, '\', '*TCs_cellpose.mat'])) % proceed if multisess cellpose time course not exist
        disp('sess cellpose time course exists, skip to next set:')
        disp(i+1)
        continue
    end

    file_name = [file_list(i).folder, '\', file_list(i).name];
    load(file_name, 'out')
    arg_ImgFolder = ['00', num2str(dataset_table(i,:).num(1))]

    [data, ~, ~, ~, ~, run_str_sess] = load_sbx_data(arg_mouse, arg_date, arg_ImgFolder);
    [outs, data_reg] = stackRegister_MA_LL(double(data), [], [], out); % re-register to get data_reg back
    clear data
    
    cd(dir_final_tif)
    tif_name = [dir_final_tif, '\cellpose_mask.mat']
    LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan'; % TODO: spare generating LL_base from function load_sbx_data
    npSub_tc = get_cellpose_timecourse(data_reg, tif_name, LL_base, arg_date, imouse, run_str_sess);
    clear data_reg
end

%%

clear global % suspect weird trace is bc residual global var affecting sbxread
clearvars -except dataset_meta nset iset dataset_table stim_type

