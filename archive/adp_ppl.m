%% adaptation analysis pipeline
%{ 
attempt to segment spagetti code into functions:
    run thru adp_plot, figure out what is needed
    generate those w functions. reuse as much as possible
    edit adp_ppl (pipeline) accordingly
%}

%% init

close all
clear
clc
cd C:\Users\lan\Documents\repos\inter\mat

%% read pre-processed mat of each dataset

dataset_list = struct;
dataset_list.mouse = [1322,1322,1322, 1323,1323,1323, 1324,1324]; % i1324 200730 LI was given up
dataset_list.date = [200803, 200804, 200806, ...
                    200720, 200721, 200723, ...
                    200728, 200729];
dataset_list.area = {'V1','LM','LI', 'V1','LM','LI', 'V1','LM'};
for iset = 1 : length(dataset_list.date)
    areamouse_seq{iset} = [dataset_list.area{1,iset} '_' num2str(dataset_list.mouse(iset))];
end

nset = length(dataset_list.date);
result_folder = cell(nset, 1);
set = cell(nset, 5);
for iset = 1 : nset
    date = num2str(dataset_list.date(iset))
    mouse = num2str(dataset_list.mouse(iset))
    imouse = ['i', mouse];
    area = dataset_list.area{1,iset}
    areamousedate = [area '_' imouse '_' date];

    result_prefix = 'C:\Users\lan\Documents\repos\inter\mat\';
    result_folder{iset} = fullfile(result_prefix, areamousedate);
    cd(result_folder{iset});

    mat_files = dir('*.mat'); 
    for imat = 1:length(mat_files) 
        set{iset, imat} = load(mat_files(imat).name); 
    end
end