%% visliz F over time: photobleaching / dying cell

% clear 
% clear all global % sbxread need to first clear global var
% close all
% clc

cd('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\002')
names = dir('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\002\*_000_000.mat');
imgMatFile = {names.name}
load(imgMatFile{1});
nframes = info.config.frames;
names = dir('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\002\*_000_000.sbx');
sbxFile = {names.name}
data_temp = sbxread(sbxFile{1}(1:end-4), 0, nframes);

cd('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\003')
names = dir('Z:\All_Staff\home\lan\Data\2P_images\i1350\220225\003\*_000_000.sbx');
sbxFile = {names.name}
data_temp2 = sbxread(sbxFile{1}(1:end-4), 0, nframes);

data = cat(4, data_temp, data_temp2);
data = squeeze(data);
data_avg = squeeze(mean(mean(data,2),1));

data_smooth = movmean(data_avg, 500);
plot(data_smooth)

cd('Z:\All_Staff\home\lan\Analysis\2P\220225_i1350')
save overall_fluorescence.mat data_avg
saveas(gcf, 'overall_fluorescence', 'jpg')
close