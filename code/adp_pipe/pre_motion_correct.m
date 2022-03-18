%%
clear all
clc

%%
database_path = 'Z:\All_Staff\home\lan\Data\2P_images\';
master_xls = [database_path, 'mat_inter\adp_dataset_master.xlsx'];
dataset_meta = readtable(master_xls);
dataset_now = dataset_meta(ismember(dataset_meta.paradigm, ...
    'bunnytop high res high lum-contrast'),:);
% dataset_now = dataset_now(dataset_now.mouse == 1369, :);
dataset_now

iset = 1
mouse = dataset_now.mouse(iset)
imouse = ['i', num2str(mouse)];
date = num2str(dataset_now.date(iset))
ImgFolder = dataset_now.num(iset); %ImgFolder = ImgFolder{1}

%%
xls_dir = fullfile(database_path, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx'); clear dataset_meta
dataset_run = readtable(xls_file.name); 
idx = find(all(ismember(dataset_run.(1),[ImgFolder,'_000_000']),2));
time = num2str(dataset_run.(8)(idx));
frame_rate = 30;

%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;


for irun = 1:nrun
    LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
    CD = [LL_base '\Data\2P_images\' imouse '\' date '\' ImgFolder(irun,:)];
%     CD = [LL_base '\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' imouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = max([temp(irun).counterValues{end}(end) info.config.frames]);

    
    fprintf(['Reading run ' num2str(irun) ', consisting of ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
%     data_temp = sbxread(imgMatFile(1,1:11),0,100000);
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cLeverUp') 
        if irun>1
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                %temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                temp(irun).cLeverDown{itrial} = temp(irun).cLeverDown{itrial}+offset;
                temp(irun).cFirstStim{itrial} = temp(irun).cFirstStim{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                if ~isempty(temp(irun).cLeverUp{itrial})
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial}+offset;
                else
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial};
                end
                if ~isempty(temp(irun).cTargetOn{itrial})
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial}+offset;
                else
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial};
                end
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
% input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc

% %% Choose register interval
nep = floor(size(data,3)./10000);
[n, n2] = subplotn(nep);
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
