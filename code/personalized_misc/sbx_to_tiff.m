%% clean up

clearvars
clear global
clc
disp('clean up done')

%% set dir

caiman_path = 'C:\Users\ll357\Documents\CaImAn\';
root_path = 'C:\Users\ll357\Documents\inter\';
% caiman_path = 'D:\Lan_temp\caiman_test_data\';
% root_path = 'C:\Users\lan\Documents\repos\inter\';

% mode = 'local' % read sbx and write tif locally
mode = 'remote' % read sbx and write tif remotely on isilon

database_path = 'Z:\All_Staff\home\lan\Data\2P_images\';
master_xls = [root_path, 'mat\adp_dataset_master.xlsx'];
dataset_meta = readtable(master_xls);
dataset_todo = dataset_meta(ismember(dataset_meta.caiman, 'todo'),:);
nset = size(dataset_todo); nset = nset(1)

%% for each sbx
time_seq = [];
for iset = 0 % test bunny 500 gcamp6s V1: 1339_210922
% 3 % test remote w i1324 V1
% 1:nset on hubel % nset:-1:1 on nuke

tic
% disp('working on dataset #')
% iset
% date = num2str(dataset_todo.date(iset))
% mouse = num2str(dataset_todo.mouse(iset)); imouse = ['i', mouse]
% area = dataset_todo.area{iset, 1}
% num = dataset_todo.num{iset, 1}

disp('convert bunny 500 gcamp6s V1')
date = '210922'
mouse = '1339'; imouse = ['i', mouse]
area = 'V1'
num = '002' % todo: 003 and 004 concat to long tif = 240K frames total

disp('prep done')

%% copy to local

cd([database_path, imouse, '\', date, '\', num, '\'])
sbx_file = [num, '_000_000.sbx'];
mat_file = [num, '_000_000.mat'];
caiman_data_path = [caiman_path, 'demos\temp_data\'];

if contains(mode, 'local')
    local_iset = [caiman_data_path, imouse, '_', date, '_', num, '\'];
    if exist(local_iset, 'dir')
        continue % if local dir exist, assume this sbx has been copied
    end
    mkdir(local_iset)
    copyfile(mat_file, local_iset);
    copyfile(sbx_file, local_iset);
    disp('copied sbx to local')
    cd(local_iset)
else 
    disp('working on sbx remotely')
end

%% load sbx

load(mat_file);
nframes = info.config.frames
data_temp = sbxread(mat_file(1,1:11), 0, nframes);
sz = size(data_temp)
disp('read sbx done')

data = permute(data_temp, [3,2,4,1]); % flip to make nrow > ncol. for easy visualiz
disp('permutation done')
data = squeeze(data);
toc

%% convert sbx to tif

% tic
disp('start saving tiff')
datetime('now')
tif_file = [imouse, '_', date, '_', num, '_multipage_100k_local.tif'];
if exist(tif_file, 'file')
    continue % if tif exist, assume this sbx has been converted
end

% saveastiff_LL(data, tif_file); 
% disp('save tiff done')
% t = toc
% time_seq(end+1) = t

fTIF = Fast_BigTiff_Write(tif_file,1,0);
tic
msg = 0;
N = nframes % /100; % tested w 1000 frames, worked. visualized in caiman
B = numel(data)*2;
for ct = 1:N
    fprintf(1,repmat('\b',[1,msg]));
    msg = fprintf(1,'%.0f/%.0f',ct,N);
    fTIF.WriteIMG(data(:,:,ct));
end
fTIF.close;
t=toc
fprintf(1,repmat('\b',[1,msg]));
fprintf(1,'\nWrite %.0f bytes in %.0f mins \n',B*N,t/60);
fprintf(1,'Write speed: %.0f MB/s \n',(B*N)/(2^20*t));

% using mapped drive isilon:
% takes 20h to convert 100K frame sbx
% takes <10h to convert 70K frame sbx? with a tifflib error in the middle ("Unable to write the current directory.")
% takes 2.5h to convert 35K frame sbx

% using local drive to read and write:
% takes 77h to convert 100K frame sbx??? why??? inspect saveastiff func
% takes 53h to convert 100K frame sbx. verified w another file. why???

%% remove sbx if local

if contains(mode, 'local')
    delete(mat_file)
    delete(sbx_file)
end
end