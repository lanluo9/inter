% testing with bunny 500 rand gcamp6s V1 dataset
% todo: generaliz

%% clean up

clear all
clear all global
clc
tic

%% set dir

caiman_path = 'C:\Users\ll357\Documents\CaImAn\';
root_path = 'C:\Users\ll357\Documents\inter\';
database_path = 'Z:\All_Staff\home\lan\Data\2P_images\';

% master_xls = [root_path, 'mat\adp_dataset_master.xlsx'];
% dataset_meta = readtable(master_xls);
% dataset_todo = dataset_meta(ismember(dataset_meta.caiman, 'todo'),:);
% dataset_todo = dataset_todo(4:end,:); % waiting for AWS
% nset = size(dataset_todo); nset = nset(1)
% 
% for iset = 1:nset

disp('working on dataset #')
% iset
% date = num2str(dataset_todo.date(iset))
% mouse = num2str(dataset_todo.mouse(iset)); imouse = ['i', mouse]
% area = dataset_todo.area{1,iset}
% num = dataset_todo.num{1,iset}

date = '210922'
mouse = '1339'; imouse = ['i', mouse]
area = 'V1'
num = '002' % todo: 003 and 004 concat to long tif = 240K frames total

%% copy to local

cd([database_path, imouse, '\', date, '\', num, '\'])
sbx_file = [num, '_000_000.sbx'];
mat_file = [num, '_000_000.mat'];
caiman_data_path = [caiman_path, 'demos\temp_data\'];
local_iset = [caiman_data_path, imouse, '_', date, '_', num, '\'];
mkdir(local_iset)
copyfile(mat_file, local_iset);
copyfile(sbx_file, local_iset);
disp('copied sbx to local')

%% prep sbx data

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

%% test parfor

parpool(3)
parfor i=1:3, c(:,i) = eig(rand(100)); end

%% convert sbx to tif locally

tic
disp('saving tiff')
options.big = true;

n = 10; % paralleliz using n workers
parpool(n)
parfor i=1:n
    frame_range = (i-1)*nframe/n : i*nframe/n
    data_chunk = data(:,:, frame_range);
    tif_file = [imouse, '_', date, '_', num, '_', num2str(n), ...
        '_multipage_100k_local.tif'];
    saveastiff(data_chunk, tif_file, options);
end


disp('save tiff done')
toc

% using mapped drive isilon:
% takes 20h to convert 100K frame sbx
% takes <10h to convert 70K frame sbx? with a tifflib error in the middle ("Unable to write the current directory.")
% takes 2.5h to convert 35K frame sbx

% using local drive to read and write:
% takes 77h to convert 100K frame sbx??? why??? inspect saveastiff func

%% remove sbx

delete(mat_file)
delete(sbx_file)

% end