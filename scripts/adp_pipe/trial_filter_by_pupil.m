%%

close all
clear all
clc

% mouse = 'i484';
% date = '220107';
% run = '001';
% time = '1029';
% data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\camaron';
% CD = [data_pn '\Data\2P_images\' mouse '\' date '\' run];
% cd(CD);
% fn = [run '_000_000_eye.mat'];

mouse = 'i1375';
date = '220915';
% run = '001'; % try with ret9pos data
% time = '1402';
run = '003'; % mix50 session 1
time = '1421';
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
CD = [data_pn '\Data\2P_images\' mouse '\' date '\' run]; % Z:\All_Staff\home\lan\Data\2P_images\i1375\220915\003
cd(CD);
fn = [run '_000_000_eye.mat'];

%%
% load data
data_temp = load(fn);
data_temp = squeeze(data_temp.data);

% crop frames to match mworks data
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);
nFrames = input.counterValues{end}(end);
data = data_temp(:,:,1:nFrames);      % the raw images...

% data = data(:, :, 1:1000); 
% disp('using test data only')

%% Crop image to isolate pupil 
% select an area as small as possible that contains pupil, but excludes bright spots
[data_crop, rect] = cropEyeData(data); % select rectangle, hit enter if happy, z to redo

%% measure pupil position/diameter
rad_range = [3 13]; % histogram of pupil radius should show both tails, otherwise adjust rad_range
Eye_data = extractEyeData(data_crop, rad_range); % check if pupil not found in too many frames
% must use full data (not truncated test data) below, to match mworks input

%% align to stimulus presentation
[rad, centroid] = alignEyeData(Eye_data,input);

%% visliz

%% discard trial by pupil radius or position

