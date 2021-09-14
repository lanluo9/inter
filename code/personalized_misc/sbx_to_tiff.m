clear all
clc


cd('Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\')
imgMatFile = '002_000_000.mat';
load(imgMatFile);
nframes = info.config.frames;
data_temp = sbxread(imgMatFile(1,1:11), 0, nframes);

sz = size(data_temp);
% data = permute(data_temp, [2,3,4,1]);
data = permute(data_temp, [3,2,4,1]); % flip to make nrow > ncol. for easy visualiz
size(data)

data_chunk = data(:, :, 1:1200); % caiman demo data = 3000 frames
% % fixme: chunk writetif -> concat tif in python
writetiff(data_chunk, ...
    'Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\002_vertical.tif')

% imshow(squeeze(data_chunk(:,:,1)))
