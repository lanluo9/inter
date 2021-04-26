%%

% cd('D:\OneDrive - Duke University\RRR增量备份\Grad\Rotation\LG_lab\archive\GaborWavelet\grating')
% ls *.png
% ori_list = [0:22.5:157.5];

for i = 1:8
    
pic_list = {'0.png' '135.png' '22.5.png' '67.5.png' '112.5.png' '157.5.png' '45.png' '90.png' };
pathname = 'D:\OneDrive - Duke University\RRR增量备份\Grad\Rotation\LG_lab\archive\GaborWavelet\grating\';
cd(pathname)
filename = pic_list{i};
im2 = imread([pathname filename]);
disp(filename)
cd ..

if size(im2,3) == 3, im2 = rgb2gray(im2); end
if size(im2,1) ~= size(im2,2)
    disp('crop image to square')
    crop_size = min(size(im2,1), size(im2,2));
    im2 = im2(1:crop_size, 1:crop_size);
end

imageSize = [1 1] * 2^7;
%%%
% This process requires Image Processing Toolbox
if size(im2,1) ~= imageSize(1), im2 = imresize(im2,imageSize); end
%%%

im2 = fliplr(im2');
im2 = double(im2);
meanValue = mean(im2(:));
im2 = im2 - meanValue;

% imageSize = size(im2); 

%% set parameters of Gabor wavelet
% param
N  = 1;               % sampling points per octave
a0 = 2^(1/N);   % scaling factor
b0 = 1;         % the unit spatial interval

param.phai = 1.5;     % band width of gabor [octave]
param.aspectRatio = 1;% aspect ratio of the gabor filter

m = ceil(log2(imageSize(1)/2));
param.m = m*N;        % the number of scale (in 6 octaves)
K = 8;                % the number of sampling orientation
param.theta0 = pi/K;  % the step size of each angular rotation

param.N = N;
param.K = K;
param.a0 = a0;
param.b0 = b0;


% tx = 1:imageSize(1); % the width of the image
% ty = 1:imageSize(2); % the height of the image
% [x y] = meshgrid(tx,ty); x = x'; y = y';

step = 0;
steps = (m+1)*K;
h = waitbar(0, 'Now analyzing...');

% set the filter bank
for ii = 0: m
    for ll = 0: K-1
        
        filterSize = 4 * 2^ii;
        tx = 1:filterSize;
        ty = 1:filterSize;
        % tx = 1:filterSize*2;
        % ty = 1:filterSize*2;
        [x,y] = meshgrid(tx,ty); x = x'; y = y';
        ctr = ( 4 + 2^(-ii) ) / 2;
        % ctr = ( 8 + 2^(-ii) ) / 2;

        GWfilter(ii+1,ll+1).even= GaborWavelet(x,y,ii,ctr,ctr,ll,param,0);
        GWfilter(ii+1,ll+1).odd = GaborWavelet(x,y,ii,ctr,ctr,ll,param,1);
    end
end


% tic
for ii = 0: m
    
    for ll = 0: K-1
        
        
        if ii == 0
            even = myConv2(GWfilter(ii+1,ll+1).even, im2, 2^ii);
            odd  = myConv2(GWfilter(ii+1,ll+1).odd , im2, 2^ii);
        else        
            even = myConv2(GWfilter(ii+1,ll+1).even, im2, 2^ii*3/2);
            odd  = myConv2(GWfilter(ii+1,ll+1).odd , im2, 2^ii*3/2);
        end
        
        step = step + 1;
        waitbar(step / steps)
        res(ii+1, ll+1).even = even;
        res(ii+1, ll+1).odd  = odd;
    end
end
% toc
param.m = param.m+1;

close(h);

%% feature value; even+odd at same scale & ori

feav = zeros(7,8,2);
for scale = 1:7
    for ori = 1:8
%         for phase = 1:2
            feav(scale, ori, 1) = sum(res(scale, ori).even(:));
            feav(scale, ori, 2) = sum(res(scale, ori).odd(:));
%         end
    end
end
subplot(1,2,1);imagesc(squeeze(feav(:,:,1))); colorbar
subplot(1,2,2);imagesc(squeeze(feav(:,:,2))); colorbar
set(gcf, 'Position', get(0, 'Screensize'));

cd('D:\OneDrive - Duke University\RRR增量备份\Grad\Rotation\LG_lab\archive\GaborWavelet\res')
save([filename(1:end-4), '.mat'], 'res', 'feav')
saveas(gcf, [filename(1:end-4), ' feav', '.jpg'])

%%
% figure
% for scale = 1:7
%     for ori = 1:8
%         subplot(7,8, 7*(scale-1)+ori)
%         imagesc(res(scale, ori).even)
%     end
% end
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, [filename(1:end-4), ' even', '.jpg'])
% 
% figure
% for scale = 1:7
%     for ori = 1:8
%         subplot(7,8, 7*(scale-1)+ori)
%         imagesc(res(scale, ori).odd)
%     end
% end
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, [filename(1:end-4), ' odd', '.jpg'])

close all
clear

end

%%

mat_list = {'0.mat' '135.mat' '22.5.mat' '67.5.mat' '112.5.mat' '157.5.mat' '45.mat' '90.mat' };

nfea = 7*8*2; nstim = 8;
F_stim = zeros(nfea, nstim);

for i=1:8
    load(mat_list{i})
    F_stim(:,i) = feav(:);
end
% imagesc(F); colorbar
