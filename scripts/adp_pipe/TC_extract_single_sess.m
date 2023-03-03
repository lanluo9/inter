%% readme
% this script register a single-session dataset, then extract timecourse from cells segmented by cellpose
% code is inherited from batch_cellpose.m and TC_extract_multi_sess.m (previously named adp_mag_mix50.m)

clear all
clearvars –global
close all
clc

% write text output from command window to log file
cd('C:\Users\ll357\Documents\inter\data\')
diary matlab_log % batch cellpose grat 8ori 3isi data in LM and LI
% diary on

%% setup

dir_meta = 'Z:\All_Staff\home\lan\Data\2P_images\mat_inter/adp_dataset_master.xlsx';
dataset_meta = readtable(dir_meta);

stim_type = 'grating' % grat_8ori_3isi
dataset_table = dataset_meta(strcmp(dataset_meta.paradigm, stim_type), :);

area_bool = logical(strcmp(dataset_table.area, 'LM') + strcmp(dataset_table.area, 'LI'));
dataset_table = dataset_table(area_bool, :);
sum(strcmp(dataset_table.area, 'LM')) % count LM data, grat_8ori_3isi
sum(strcmp(dataset_table.area, 'LI')) % count LI

nset = size(dataset_table,1);

%% draw cellpose_tif for each data

for iset = 1:nset
% iset = 1;

iset, nset
dataset_now = dataset_table(iset,:);
arg_mouse = dataset_now.mouse
arg_date = num2str(dataset_now.date)
arg_ImgFolder = dataset_now.num{1}

ref_ImgFolder = arg_ImgFolder; % dataset_table(1,:).num{1} for multisess, ref is the same as arg_ImgFolder if single sess
run_str_ref = ['runs-' ref_ImgFolder(1,:)]

imouse = ['i', num2str(arg_mouse)];
dir_analysis = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse, ...
    '\', arg_date, '_', imouse, '_runs-', arg_ImgFolder];
try cd(dir_analysis)
catch
    mkdir(dir_analysis)
    cd(dir_analysis)
end

% check if time course exists
if ~isempty(dir('*TC*.mat'))
    disp('timecourse exists, skip to next dataset:')
    continue
end

disp('add your stim_type into function: bunny, mix, grat6, grating(grat1)')
[data_reg, ~, date, imouse, run_str] = get_data_reg_cellpose_tif(... % register every sess against itself: single sess registration
    arg_mouse, arg_date, arg_ImgFolder, stim_type, run_str_ref); 
disp(['got data_reg & cellpose tif for ', num2str(arg_mouse), ' ', arg_date, ' ', arg_ImgFolder])
close all

disp('wait for batch_cellpose.ipynb')
while ~exist('cellpose_mask.mat','file')
    pause(60) % wait for cellpose to run in ipynb
end
disp('got cellpose mask, now extract TC')

tif_name = 'cellpose_mask.mat';
LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
npSub_tc = get_cellpose_timecourse(data_reg, tif_name, ...
    LL_base, arg_date, imouse, run_str);

disp([num2str(iset), ' set done out of ', num2str(nset)])
clear data_reg
clear global % suspect weird trace is bc residual global var affecting sbxread

end