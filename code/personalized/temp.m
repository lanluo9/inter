%%

% subplot(1,2,1)
% imagesc(data_avg)
% subplot(1,2,2)
% imagesc(reg_img_mean)
% set(gcf, 'Position', get(0, 'Screensize'));

%% get path names
clear all
clc

date = '200803';
mouse = 'i1322';
ImgFolder = strvcat('002');
time = strvcat('1140');
frame_rate = 30;

doFromRef = 0;
% ref = strvcat('001'); % what is ref?
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);
LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LL_base '\Data\2P_images\' mouse '\' date '\' ImgFolder(irun,:)];
%     CD = [LL_base '\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = max([temp(irun).counterValues{end}(end) info.config.frames]);

    
    fprintf(['Reading run ' num2str(irun) ', consisting of ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
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
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc

%%
cd C:\Users\lan\Documents\repos\inter\code\py_playground
load suite2p_reg_mean_img.mat

[out, data_reg] = stackRegister(data, reg_img_mean);    % out=shift

%%
cd Z:\All_Staff\home\lan\Analysis\2P\200803_1322\200803_1322_runs-002
load 200803_i1322_runs-002_mask_cell_addfake.mat

cell_mask_orig = mask_np;
[out2, cell_mask_align_w_suite2p] = stackRegister_MA(cell_mask_orig, [],[],out);  % add shift to cell mask


[out3, cell_mask_align_w_suite2p_2] = stackRegister_MA(mask_cell, [],[],out);
imagesc(cell_mask_align_w_suite2p_2)

save mask_cell_shift_align.mat cell_mask_align_w_suite2p_2 out out2 out3

%%

tt = cell_mask_align_w_suite2p_2;
% histogram(cell_mask_align_w_suite2p_2(:),10)


tt2 = tt;
tt2(tt2<3) = nan;

tt2(~isnan(tt2)) = 1;
imagesc(tt2)

% save mask_cell_shift_flat.mat tt2