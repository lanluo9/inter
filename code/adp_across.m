%% set up directory

close all
clear
clc

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
ll_fn = fullfile(fn_base, 'home\lan'); 
data_fn = fullfile(ll_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); 
tc_fn = fullfile(ll_fn, 'Analysis\2P');

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806,...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};

%%

iset = 1;
% for iset = 1 : length(dataset_list.date)
    date = num2str(dataset_list.date(iset));
    mouse = num2str(dataset_list.mouse(iset));
% end
imouse = ['i', mouse];
xls_dir = fullfile(data_fn, imouse, date);
cd(xls_dir)
xls_file = dir('*.xlsx'); 
dataset_meta = readtable(xls_file.name); 
time = dataset_meta.(8){end};
ImgFolder = dataset_meta.(1){end}(1:3)
frame_rate = dataset_meta.(5){end};

%%
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' imouse];
datemouserun = [date '_' imouse '_' run_str];

fName = fullfile(mworks_fn, ['data-' imouse '-' date '-' time '.mat']);
load(fName); % load behavior data "input"

CD = fullfile(data_fn, imouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata "info"

tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs_addfake.mat']); % load time course including fake targ resp
% fix bug later: tc folder vs retinotopy folder in Analysis naming convention

result_folder = 'C:\Users\lan\Documents\repos\inter\code\';
cd C:\Users\lan\Documents\repos\inter\code\
