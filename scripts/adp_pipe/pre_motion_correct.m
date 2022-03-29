%{
pre register each run 002/003/004 by matlab
-> mat write to tiff
-> caiman segment -> caiman multisess
try with gcamp8f V1?
%}

%%% prep
clear all
clc

database_path = 'Z:\All_Staff\home\lan\Data\2P_images\';
master_xls = [database_path, 'mat_inter\adp_dataset_master.xlsx'];
dataset_meta = readtable(master_xls);
dataset_now = dataset_meta(ismember(dataset_meta.paradigm, ...
    'bunnytop high res high lum-contrast'),:);
% dataset_now = dataset_now(dataset_now.mouse == 1369, :);
dataset_now
nset = size(dataset_now); nset = nset(1)

for iset = 1:nset

mouse = dataset_now.mouse(iset)
imouse = ['i', num2str(mouse)];
date = num2str(dataset_now.date(iset))
ImgFolder = dataset_now.num(iset); ImgFolder = ImgFolder{1}

xls_dir = fullfile(database_path, imouse, date); cd(xls_dir)
xls_file = dir('*.xlsx'); clear dataset_meta
dataset_run = readtable(xls_file.name); 
idx = find(all(ismember(dataset_run.(1),[ImgFolder,'_000_000']),2));
time = num2str(dataset_run.(8)(idx));
frame_rate = 30;

nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);

%%% load
for irun = 1:nrun
    LL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lan';
    CD = [LL_base '\Data\2P_images\' imouse '\' date '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    
    stim_input_mat = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' imouse '-' date '-' time(irun,:) '.mat'];
    load(stim_input_mat);
    
    nframes = max([temp(irun).counterValues{end}(end) info.config.frames]);
    fprintf(['Reading run ' num2str(irun) ', consisting of ' num2str(min(nframes)) ' frames \r\n'])
    data = sbxread(imgMatFile(1,1:11),0,min(nframes));
    data = squeeze(data);
end

%%% Choose register interval
nep = floor(size(data,3)./10000);
[n, n2] = subplotn(nep);
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nep
    subplot(n,n2,i)
    imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3))
    title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))])
end

%%% select ref episode to register against

clear input
prompt = "select which episode: ";
select = input(prompt,"s");
select = str2num(select)
if isempty(select)
    select = 4;
end

start_idx = select * 10000 + 1;
stop_idx = select * 10000 + 500;
data_avg = mean(data(:,:,start_idx:stop_idx),3);

%%% register
if exist(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]))
    load(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']))
%     save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
else
    tic; [out, data_reg] = stackRegister(data,data_avg); toc;
    mkdir(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str]))
    save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
%     save(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_input.mat']), 'input')
end
size(data_reg)

%%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
set(gcf, 'Position', get(0, 'Screensize'));
print(fullfile(LL_base, 'Analysis\2P', [date '_' imouse], [date '_' imouse '_' run_str], [date '_' imouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')
close all

%%% write mat var to tiff

data_full = permute(data_reg, [2,1,3]); % flip to make nrow > ncol. for easy visualiz
disp('permutation done')
data_size = size(data_full)
nframes = data_size(3)

cd ..
disp('start saving tiff')
datetime('now')
num = ImgFolder;
tif_file = [imouse, '_', date, '_', num, '_registered.tif'];

fTIF = Fast_BigTiff_Write(tif_file,1,0);
tic
msg = 0;
N = nframes;
B = numel(data_full)*2;
for ct = 1:N
    fprintf(1,repmat('\b',[1,msg]));
    msg = fprintf(1,'%.0f/%.0f',ct,N);
    fTIF.WriteIMG(data_full(:,:,ct));
end
fTIF.close;
fprintf(1,repmat('\b',[1,msg]));
t=toc
fprintf(1,'\nWrite %.0f bytes in %.0f mins \n',B*N,t/60);
fprintf(1,'Write speed: %.0f MB/s \n',(B*N)/(2^20*t));

end