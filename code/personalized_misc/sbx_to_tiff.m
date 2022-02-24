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
% dataset_now = dataset_meta(ismember(dataset_meta.caiman, 'todo'),:);
dataset_now = dataset_meta(ismember(dataset_meta.paradigm, 'bunnytop'),:);
dataset_now = dataset_now(ismember(dataset_now.area, 'LM'),:);
dataset_now
nset = size(dataset_now); nset = nset(1)

% data_seq = 'single';
% if sum(ismember(dataset_now.paradigm, 'bunny500'))
%     data_seq = 'sequential';
% end

%% for each sbx

% if contains(data_seq, 'sequential')
%     data_full = [];
% end

for iset = 1:nset
    
tic
disp('working on dataset #')
iset
date = num2str(dataset_now.date(iset))
mouse = num2str(dataset_now.mouse(iset)); imouse = ['i', mouse]
area = dataset_now.area{iset, 1}
num = dataset_now.num{iset, 1}

disp('prep done')

%% copy to local or work remotely

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
% data = data(:,:,1:70000); % get first 70k frames = 1 round of bunny500

% data_full = cat(3,data_full,data); % concat 3 recordings of bunny 500
data_full = data;
data_size = size(data_full)
nframes = data_size(3)
toc

% %% save max proj: pre-registration movie
% 
% max_proj = max(data,[],3);
% size(max_proj)
% 
% imagesc(max_proj)
% set(gcf, 'Position',  [0, 0, width(max_proj), height(max_proj)])

%% convert sbx to tif

% tic
cd ..
disp('start saving tiff')
datetime('now')
tif_file = [imouse, '_', date, '_', num, '_multipage.tif'];
% if exist(tif_file, 'file')
%     continue % if tif exist, assume this sbx has been converted
% end

fTIF = Fast_BigTiff_Write(tif_file,1,0);
tic
msg = 0;
N = nframes;
B = numel(data_full)*2;
for ct = 1:N
    fprintf(1,repmat('\b',[1,msg]));
    msg = fprintf(1,'%.0f/%.0f',ct,N);
    fTIF.WriteIMG(data_full(:,:,ct));
end
fTIF.close;
fprintf(1,repmat('\b',[1,msg]));
t=toc
fprintf(1,'\nWrite %.0f bytes in %.0f mins \n',B*N,t/60);
fprintf(1,'Write speed: %.0f MB/s \n',(B*N)/(2^20*t));

%% remove sbx if local

if contains(mode, 'local')
    delete(mat_file)
    delete(sbx_file)
end

end