%% 

clear all
clc

%% 

database_path = 'C:/Users/GlickfeldLab/Documents/test/inter/';
master_xls = [database_path, 'data/batch_cellpose.csv'];
dataset_meta = readtable(master_xls);
dataset_bunny = dataset_meta(ismember(dataset_meta.stim_type, 'bunny'),:);
dataset_grat = dataset_meta(ismember(dataset_meta.stim_type, 'grating'),:);

%%

nset = size(dataset_bunny,1);
disp('cellpose segment bunny datasets first')

for iset = 1 : nset

iset, nset
dataset_now = dataset_bunny(iset,:);
arg_mouse = dataset_now.mouse
arg_date = num2str(dataset_now.date)
arg_ImgFolder = dataset_now.num(1); arg_ImgFolder = ['00', num2str(arg_ImgFolder)]
disp('for bunny data, take sess 002 only to avoid repetitive depth')

% check if cellpose time course exists
imouse = ['i', num2str(arg_mouse)];
dir_analysis = ['Z:\All_Staff\home\lan\Analysis\2P\', arg_date, '_', imouse, ...
    '\', arg_date, '_', imouse, '_runs-', arg_ImgFolder];
try cd(dir_analysis)
catch
    mkdir(dir_analysis)
    cd(dir_analysis)
end
if ~isempty(dir('*TCs_cellpose.mat'))
    disp('cellpose time course exists, skip to next set:')
    disp(iset+1)
    continue
end

[data_reg, LL_base, date, imouse, run_str] = get_data_reg_cellpose_tif(...
    arg_mouse, arg_date, arg_ImgFolder);
disp('got data_reg, waiting for cellpose')

tic
while ~exist('cellpose_mask.mat','file')
    % pwd should be like Z:\All_Staff\home\lan\Analysis\2P\220310_i1369\220310_i1369_runs-002
    pause(60) % wait for cellpose to run in ipynb
    toc
end

disp('got cellpose mask, extracting TC')
npSub_tc = get_cellpose_timecourse(data_reg, ...
    LL_base, arg_date, imouse, run_str);

disp([num2str(iset), ' set done out of ', num2str(nset)])
end