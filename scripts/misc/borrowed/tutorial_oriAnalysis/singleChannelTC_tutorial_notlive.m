%% prepare folder names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey'); 
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); % mwork = behavior data
fnout = fullfile(lg_fn, 'Analysis\2P\test'); % output in analysis folder

%% specify file names & load data
date = '200118';
ImgFolder = '002';
time = strvcat('1508'); % catenate multiple time vertically to get a matrix of file/folder name
mouse = 'i1312';
frame_rate = 15.5;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName); % load behavior data, aka "input"
CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata, aka "info" % check content

%% visualize episodes
totframes = input.counterValues{end}(end); % total # of frames
fprintf(['Reading ' num2str(totframes) ' frames \r\n'])
tic; data = sbxread([ImgFolder '_000_000'], 0, totframes); disp('finished reading'); toc;
fprintf(['Data dimension: ' num2str(size(data)) '\n'])
data = squeeze(data);
fprintf(['Data new dimension: ' num2str(size(data)) '\n'])
nframes = 500; % nframes to average
nskip = 1500; % nframes to skip for each average. avg 500 frame -> skip 1k frame -> etc

% choose target both stable and at middle of stack, to account for x-y shift

%%
nep = floor(size(data,3)./nskip); % # of subplot. ep = episode
[n, n2] = subplotn(nep); 
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep; 
    subplot(n,n2,i); 
    imagesc(mean(data(:,:,1+((i-1)*nskip):nframes+((i-1)*nskip)), 3)); % take avg over nframes
%     title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip))]); 
    title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip)) '-' num2str(i)]); 
end

f=gcf;
w = waitforbuttonpress; %click on subplot
if w == 0 % Returns 0 when terminated by mouse, 1 when keypress
    axesClicked = gca;
    allAxes = findobj(f.Children,'Type','axes');
    numClicked = find(axesClicked==allAxes);
    close all
end
fprintf(['Selected subplot ' num2str(numClicked) '\n']) % i pressed subplot 6 and numClicked = 4. why? -> test this
% which one should i click? cells look sharpest. could plot a smaller zoomed-in fov to make it more obvious
% if vasculature unfocus, there is drift in Z which cannot be rectified

%% analyze & save specific episode
data_avg = mean(data(:,:,1+((numClicked-1)*nskip):nframes+((numClicked-1)*nskip)),3); % avg of the clicked 500-frame episode
[out, data_reg] = stackRegister(data, data_avg); % data as stack, data_avg as target
data_reg_avg = mean(data_reg,3); % avg over all frames of registered data

mkdir(fullfile(fnout, datemouse, datemouserun)) % make analysis result folder
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg') 
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input') % input = behavior data

ind = [1 nep];
for i = 1:length(ind) 
    subplot(2,1,i); 
    ix = ind(i);
    imagesc(mean(data_reg(:,:,1+((ix-1)*nskip):nframes+((ix-1)*nskip)),3)); 
    title([num2str(1+((ix-1)*nskip)) '-' num2str(nframes+((ix-1)*nskip))]); 
end
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_first&last.pdf']), '-dpdf')
close
imagesq(data_reg_avg); 
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']), '-dpdf') % FOV = Field Of View

%% normalize & sort by dir & filter
clear data
% trial = no-stim (nScansOff) period + stim present (nScansOn)
nOn = input.nScansOn; % behavior data "input"
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2); %this is a cell array with one grating direction value per trial, so length = ntrials
sz = size(data_reg);
data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]); % data ordered by trial
fprintf(['Size of data_tr is ' num2str(size(data_tr))])
data_f = mean(data_tr(:,:,nOff/2:nOff,:), 3); % (y,x,frame,dir). nOff/2 account for ca reaction time
data_df = bsxfun(@minus, double(data_tr), data_f); % subtract avg (baseline) of trial
data_dfof = bsxfun(@rdivide, data_df, data_f); % normalize by baseline. dfof = delta_f/f0 = (f-f0)/f0

clear data_f data_df data_tr
Dir = celleqel2mat_padded(input.tGratingDirectionDeg); %transforms cell array into matrix (1 x ntrials)
Dirs = unique(Dir);
nDirs = length(Dirs);
data_dfof_avg = zeros(sz(1),sz(2),nDirs); % avg over grating directions
%create empty matrix with FOV for each direction: nYpix x nXPix x nDir

nStim = nDirs;
[n, n2] = subplotn(nDirs); %function to optimize subplot number/dimensions
figure('units','normalized','outerposition',[0 0 1 1]);
for idir = 1:nDirs
    ind = find(Dir == Dirs(idir)); %find all trials with each direction
    data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4); %average over On frames and then over same-dir trials
    subplot(n,n2,idir)
    imagesc(data_dfof_avg(:,:,idir))
end
clear data_dfof

myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = imfilter(data_dfof_avg, myfilter); % dfof by dir, filtered by lowpass gaussian
data_dfof_max = max(data_dfof_avg_all,[],3); %finds all active cells by taking max projection across dirs
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')

%% get time course of neuropil-adjusted cell 
data_dfof = cat(3,data_dfof_max, data_dfof_avg_all);
mask_data = data_dfof;
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0; %blacks out old cells
    bwout = imCellEditInteractiveLG_LL(mask_data_temp); %selection GUI
    mask_all = mask_all + bwout; %adds new cells to old cells
    mask_exp = imCellBuffer(mask_all, 3) + mask_all; %creates buffer around cells to avoid fusing
    close all
end % select all bright cells
mask_cell = bwlabel(mask_all); %turns logical into numbered cells

figure;
imagesc(mask_cell)
t = string(datetime('now'));
t = replace(t, ':', '_'); t = replace(t, ' ', '_');
saveas(gcf, ['mask_cell_', t, '.png'])


nMaskPix = 5; %thickness of neuropil ring in pixels -> neighbor's contamination
nBuffPix = 3; %thickness of buffer between cell and ring -> account for illumination of surrounding of the cell itself
mask_np = imCellNeuropil(mask_cell, nBuffPix, nMaskPix); % account for no-cell area around cell
% point spread function -> contamination from neighbor cell
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
% clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

data_tc = stackGetTimeCourses(data_reg, mask_cell); %apply cell mask to img stack to get time courses
fprintf(['data_tc is ' num2str(size(data_tc))]) 
nCells = size(data_tc,2);
down = 5; %number of frames to average
data_reg_down = stackGroupProject(data_reg,down); %averages every 5 frames in stack  
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell); % downsampled time course

sz = size(data_reg);
np_tc = zeros(sz(3), nCells); % np_tc = time course of neuropil-adjusted cell mask
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg, mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down, mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '\n']) 
end

%% minimize neuropil contribution
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down - tcRemoveDC(np_tc_down*ii(i))); % tcRemoveDC = tc - mean (direct current)
    % assume sparse ca signal -> long tail & positive skewness
end
[max_skew, ind] =  max(x,[],1); % maximize skew by i value
np_w = 0.01*ind; % np_weight
npSub_tc = data_tc - bsxfun(@times, tcRemoveDC(np_tc), np_w);
clear data_reg data_reg_down
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

