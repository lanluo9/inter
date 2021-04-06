function [npSub_tc, frame_rate, input_behav, info, result_folder] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area)
% function [frame_rate, input_behav, info, result_folder] = load_xls_tc_stim(data_fn, mworks_fn, tc_fn, date, imouse, area)

% input: 
%     "data_fn" lab note xls
%     "mworks_fn" stim input
%     "tc_fn" of 2P img timecourse
%      date, imouse, area
% 
% output:
%     result_folder directory
%     input_behav for stim input
%     info for 2P img metadata
%     npSub_tc for timecourse nframe x ncell

xls_dir = fullfile(data_fn, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx'); clear dataset_meta
dataset_meta = readtable(xls_file.name); 
time = dataset_meta.(8)(end);
ImgFolder = dataset_meta.(1){end}(1:3);
frame_rate = dataset_meta.(5)(end);

run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' imouse]; datemouserun = [date '_' imouse '_' run_str];
areamousedate = [area '_' imouse '_' date];

fName = fullfile(mworks_fn, ['data-' imouse '-' date '-' num2str(time) '.mat']);
temp = load(fName); % load behavior data "input"
input_behav = temp.input; clear temp

CD = fullfile(data_fn, imouse, date, ImgFolder); cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata "info"
tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs_addfake.mat']); 
% fix bug later: retinotopy folder in Analysis naming convention should adhere to tc folder

result_prefix = 'C:\Users\lan\Documents\repos\inter\mat\';
result_folder = fullfile(result_prefix, areamousedate); if ~exist(result_folder); mkdir(result_folder); end
cd(result_folder);
