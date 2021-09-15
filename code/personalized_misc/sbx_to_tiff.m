clear all
clc

cd('Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\')
imgMatFile = '002_000_000.mat';
load(imgMatFile);
nframes = info.config.frames;
data_temp = sbxread(imgMatFile(1,1:11), 0, nframes);

sz = size(data_temp);
data = permute(data_temp, [3,2,4,1]); % flip to make nrow > ncol. for easy visualiz

% to write large tiff, matlab - pref - general - java heap memory (make larger)
writetiff(data, 'Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\002_vertical_full.tif')
% writetiff for 100K frames takes about 1 hour
