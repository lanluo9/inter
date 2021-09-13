% open sbxread.m
% open writetiff.m

cd('Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\')
% fname = '002_000_000.sbx';
% sbx2tif(fname)

imgMatFile = '002_000_000.mat';
load(imgMatFile);
nframes = info.config.frames;
data_temp = sbxread(imgMatFile(1,1:11), 0, nframes);

sz = size(data_temp);
% sz = [sz(2:4), sz(1)];
% data = reshape(data_temp, sz);
data = permute(data_temp, [2,3,4,1]);
size(data)

data_chunk = data(:, :, 1:1000);
writetiff(data_chunk, 'Z:\All_Staff\home\lan\Data\2P_images\i1329\201209\002\002.tiff')