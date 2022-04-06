t = load('Z:\All_Staff\home\lan\Analysis\2P\220310_i1369\220310_i1369_runs-002\data_dfof.mat');
stim_resp_gauss = t.data_dfof_max_gauss;
stim_resp_gauss = stim_resp_gauss';

tif_file = ['stim_resp_gauss.tif'];

fTIF = Fast_BigTiff_Write(tif_file,1,0);
tic
msg = 0;
N = 1;
B = numel(stim_resp_gauss)*2;
for ct = 1:N
    fprintf(1,repmat('\b',[1,msg]));
    msg = fprintf(1,'%.0f/%.0f',ct,N);
    fTIF.WriteIMG(stim_resp_gauss(:,:,ct));
end
fTIF.close;
fprintf(1,repmat('\b',[1,msg]));
t=toc
fprintf(1,'\nWrite %.0f bytes in %.0f mins \n',B*N,t/60);
fprintf(1,'Write speed: %.0f MB/s \n',(B*N)/(2^20*t));

%%

% unique(mask_cell(:))
% unique(mask_np(:))

t = load('C:\Users\GlickfeldLab\Documents\test\inter\scripts\ipynb\cellpose_mask.mat');
cellpose_mask = t.cellpose_mask;
mask_cell = cellpose_mask;

%%

mask_np = imCellNeuropil(mask_cell, 3, 5);
% save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_mask_cell_addfake.mat']), 'data_dfof', 'mask_cell', 'mask_np')
% clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5; % down sampling
sz = size(data_reg);
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2)

np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s /n']) 
     disp(' ')
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew, ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], ...
    [date '_' imouse '_' run_str], ...
    [date '_' imouse '_' run_str '_TCs_cellpose.mat']), ...
    'data_tc', 'np_tc', 'npSub_tc')
clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
disp('TC extraction complete')