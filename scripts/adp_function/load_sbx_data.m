function [data, behav_input, LL_base, date, imouse, run_str] = load_sbx_data(arg_mouse, arg_date, arg_ImgFolder)

disp('start running get_data_reg_cellpose_tif.m')
mouse = arg_mouse;
imouse = ['i', num2str(mouse)];
date = num2str(arg_date);
ImgFolder = arg_ImgFolder;

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images';
xls_dir = fullfile(fn_base, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx'); clear dataset_meta
dataset_meta = readtable(xls_file.name); 
idx = find(all(ismember(dataset_meta.(1),[ImgFolder,'_000_000']),2));
time = num2str(dataset_meta.(8)(idx));
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);

% load and register
disp('start loading sbx')
tic
data = [];
clear temp
trial_n = [];
offset = 0;

for irun = 1 : nrun
    LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
    CD = [LL_base '\Data\2P_images\' imouse '\' date '\' ImgFolder(irun,:)];
    cd(CD);
    
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' imouse '-' date '-' time(irun,:) '.mat'];
    behav_mat = load(fName);
    temp(irun) = behav_mat.input;
    behav_input = behav_mat.input;
    nframes = max([temp(irun).counterValues{end}(end) info.config.frames]);
    fprintf(['Reading run ' num2str(irun) ', consisting of ' num2str(min(nframes)) ' frames \r\n'])

    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    offset = offset+min(nframes);
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
% behav_input = concatenateDataBlocks(temp);
clear data_temp temp
toc

end