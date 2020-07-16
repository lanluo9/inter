%% 

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey'); 
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data'); % mwork = behavior data
tc_fn = fullfile(lg_fn, 'Analysis\2P');
fnout = fullfile(lg_fn, 'Analysis\2P\test'); % output in analysis folder

date = '200118';
ImgFolder = '002';
time = strvcat('1508'); % catenate multiple time vertically to get a matrix of file/folder name
mouse = 'i1312';
frame_rate = 15.5;
run_str = catRunName(ImgFolder, 1);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

%% load data

fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time '.mat']);
load(fName); % load behavior data, aka "input"

CD = fullfile(data_fn, mouse, date, ImgFolder);
cd(CD);
imgMatFile = [ImgFolder '_000_000.mat'];
load(imgMatFile); % load 2P img metadata, aka "info" % check content

tc_name = fullfile(tc_fn, datemouse, datemouserun);
load([tc_name, '\', datemouserun, '_TCs.mat']);

%% use npSub_tc to conduct further analysis

nOn = input.nScansOn; % behavior data "input"
nOff = input.nScansOff;
trial_len = nOn + nOff;
ntrials = size(input.tGratingDirectionDeg,2);
size(npSub_tc) % nframe * ncell
ncell = size(npSub_tc, 2);

tc_trials = reshape(npSub_tc,[nOn+nOff, ntrials, size(npSub_tc,2)]); 
tc_trials = permute(tc_trials, [3,2,1]);
size(tc_trials) % ncell * ntrial * trial_len

Dir = celleqel2mat_padded(input.tGratingDirectionDeg); 
convert_idx = Dir>=180;
Ori = Dir;
Ori(convert_idx) = Ori(convert_idx) - 180; 
Ori_list = unique(Ori);
nOri = length(Ori_list);

% tt = mean(mean(tc_trials,1),2);
% plot(1:trial_len, squeeze(tt)) 
% trial = 0-60 off + 61-90 on. signal decays >2/3 after 0-30 off
% as short as 3 frames would suffice. now take 10 frames as window len:

for iori = 1 : nOri
    idx = find(Ori == Ori_list(iori)); 
    for icell = 1 : ncell
        base_win = 
        resp_win = 
        
    end
end







