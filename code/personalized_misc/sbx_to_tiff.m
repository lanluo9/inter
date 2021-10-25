clear all
clear all global
clc

tic
% cd('Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\')
cd('C:\Users\ll357\Documents\CaImAn\demos\demo_data_002') % copy and convert locally!

imgMatFile = '002_000_000.mat';
load(imgMatFile);
nframes = info.config.frames;
nframes
data_temp = sbxread(imgMatFile(1,1:11), 0, nframes);
size(data_temp)
disp('read sbx done')

sz = size(data_temp);
data = permute(data_temp, [3,2,4,1]); % flip to make nrow > ncol. for easy visualiz
disp('permutation done')

% % to write large tiff, matlab - pref - general - java heap memory (make larger)
% writetiff(data, ...
%     'Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\002_vertical_full.tif')
% % writetiff for 100K frames takes about 1 hour
% % but it outputs single-page tif, cannot feed to caiman

% elisabeth:
% data = squeeze(sbxread(fname, 0, 5000));
% saveastiff(data, ‘fname.tif’);

data = squeeze(data);
% data_chop = data(:,:,1:70000);
toc

tic
disp('start saving tiff')
saveastiff(data, '002_multipage_100k_local.tif');
%     'Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\002_multipage_100k_local.tif');
% using mapped drive isilon:
% takes 20h to convert 100K frame sbx
% takes <10h to convert 70K frame sbx? with a tifflib error in the middle ("Unable to write the current directory.")
% takes 2.5h to convert 35K frame sbx
% using local drive to read and write:
% takes 77h to convert 100K frame sbx??? why???
toc