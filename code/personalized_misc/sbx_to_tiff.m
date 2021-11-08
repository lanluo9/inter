%% clean up

clearvars
clear global
clc
tic

%% set dir

caiman_path = 'C:\Users\ll357\Documents\CaImAn\';
root_path = 'C:\Users\ll357\Documents\inter\';
database_path = 'Z:\All_Staff\home\lan\Data\2P_images\';
master_xls = [root_path, 'mat\adp_dataset_master.xlsx'];
dataset_meta = readtable(master_xls);
dataset_todo = dataset_meta(ismember(dataset_meta.caiman, 'todo'),:);
% dataset_todo = dataset_todo(4:end,:); % waiting for AWS
nset = size(dataset_todo); nset = nset(1)

time_seq = [];
for iset = 1 %:nset

disp('working on dataset #')
iset
date = num2str(dataset_todo.date(iset))
mouse = num2str(dataset_todo.mouse(iset)); imouse = ['i', mouse]
area = dataset_todo.area{iset, 1}
num = dataset_todo.num{iset, 1}

%% copy to local

cd([database_path, imouse, '\', date, '\', num, '\'])
sbx_file = [num, '_000_000.sbx'];
mat_file = [num, '_000_000.mat'];
caiman_data_path = [caiman_path, 'demos\temp_data\'];
local_iset = [caiman_data_path, imouse, '_', date, '_', num, '\'];
if ~exist(local_iset, 'dir')
    continue % if local dir exist, assume this sbx has been copied
end
mkdir(local_iset)
copyfile(mat_file, local_iset);
copyfile(sbx_file, local_iset);
disp('copied sbx to local')

%% convert sbx to tif locally

cd(local_iset)
load(mat_file);
nframes = info.config.frames
data_temp = sbxread(mat_file(1,1:11), 0, nframes);
sz = size(data_temp)
disp('read sbx done')

data = permute(data_temp, [3,2,4,1]); % flip to make nrow > ncol. for easy visualiz
disp('permutation done')
data = squeeze(data);
toc

tic
disp('saving tiff')
tif_file = [imouse, '_', date, '_', num, '_multipage_100k_local.tif'];
if ~exist(tif_file, 'file')
    continue % if tif exist, assume this sbx has been converted
end
saveastiff(data, tif_file);
disp('save tiff done')
t = toc
time_seq(end+1) = t

% using mapped drive isilon:
% takes 20h to convert 100K frame sbx
% takes <10h to convert 70K frame sbx? with a tifflib error in the middle ("Unable to write the current directory.")
% takes 2.5h to convert 35K frame sbx

% using local drive to read and write:
% takes 77h to convert 100K frame sbx??? why??? inspect saveastiff func
% takes 53h to convert 100K frame sbx. why???

%% remove sbx

delete(mat_file)
delete(sbx_file)

end