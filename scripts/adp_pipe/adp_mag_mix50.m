%% 

clear all
close all
clc

%% 

% database_path = 'C:/Users/GlickfeldLab/Documents/test/inter/';
database_path = 'C:\Users\ll357\Documents\inter\';
master_xls = [database_path, 'data/mix50_grat1.csv'];
dataset_meta = readtable(master_xls);

%%

clearvars -except dataset_meta
clearvars –global

stim_type = 'mix' % 3 sessions of mix50 must be registered to align

dataset_mix = dataset_meta(ismember(dataset_meta.stim_type, 'mix'),:);
dataset_grat = dataset_meta(ismember(dataset_meta.stim_type, 'grating'),:);

if strcmp(stim_type, 'mix')
    dataset_table = dataset_mix;
elseif strcmp(stim_type, 'grating')
    dataset_table = dataset_grat
end
nset = size(dataset_table,1);
% disp('cellpose segment grating datasets first')

% for iset = 1:nset
iset = 1

iset, nset
dataset_now = dataset_table(iset,:);
arg_mouse = dataset_now.mouse
arg_date = num2str(dataset_now.date)
arg_ImgFolder = dataset_now.num(1);
arg_ImgFolder = ['00', num2str(arg_ImgFolder)]
if strcmp(stim_type, 'mix')
    disp('take first session 002 of mix50 for now')
end

% check if cellpose time course exists
imouse = ['i', num2str(arg_mouse)];
dir_analysis = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse, ...
    '\', arg_date, '_', imouse, '_runs-', arg_ImgFolder];
try cd(dir_analysis)
catch
    mkdir(dir_analysis)
    cd(dir_analysis)
end
% if ~isempty(dir('*TCs_cellpose.mat'))
%     disp('cellpose time course exists, skip to next set:')
%     disp(iset+1)
%     continue
% end

[data_reg, LL_base, date, imouse, run_str] = get_data_reg_cellpose_tif(...
    arg_mouse, arg_date, arg_ImgFolder, stim_type);
disp('got data_reg, waiting for cellpose')

%%
tic
while ~exist('cellpose_mask.mat','file')
    % pwd should be like \Analysis\2P\220310_i1369\220310_i1369_runs-002
    pause(60) % wait for cellpose to run in ipynb
    toc
end

disp('got cellpose mask, extracting TC')
npSub_tc = get_cellpose_timecourse(data_reg, ...
    LL_base, arg_date, imouse, run_str);

disp([num2str(iset), ' set done out of ', num2str(nset)])
clear global
clearvars -except dataset_meta nset iset dataset_table stim_type
% suspect weird trace is bc residual global var affecting sbxread

% end

%% register mix50 multisession

data_reg_002 = data_reg;
clear data_reg
