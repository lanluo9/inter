clc

n = 500
tic
A = ([1 0 1; -1 -2 0; 0 1 -1]);
for i = 1:n
    e = eig(A);
    % disp(e)
end
toc

tic
A = gpuArray([1 0 1; -1 -2 0; 0 1 -1]);
for i = 1:n
    e = eig(A);
    % disp(e)
end
toc

%%

% delete(gcp('nocreate'))
% parpool('local',4)
% 
% parfor i = 1:100
%     a = 0
%     a = a*4 - 1000
% end

%%

clear all; close all; clc; fclose all;
%generate some data
N=1E3;
IM=imread('landOcean.jpg');
IM = uint16(sum(IM,3));
IM = IM(100:310,960:1170);
IM = IM-min(IM(:));
IM=IM*(2^15/max(IM(:)));
IM = repmat(IM,[1,1,N])+randi((2^15)-1,[size(IM,1),size(IM,2),N],'uint16');
S = (numel(IM)/N*2)/2^20;


%imread writespeed
methods = {'imwrite','tifflib','fTIF'};
for M = 1:length(methods)
    method = methods{M};
    %file
    filename = [method,'.tif'];
    if exist(filename,'file'), delete(filename);end
    switch method
        case 'imwrite'
            %timing vector
            t = zeros(1,100+1);
            tic;
            imwrite(IM(:,:,1),filename);
            t(2)=toc;
            for ct = 2:100
                imwrite(IM(:,:,ct),filename,'WriteMode','append');
                t(ct+1)=toc;
            end
        case 'tifflib'
            %timing vector
            t = zeros(1,200+1);
            tic;
            tf = Tiff(filename,'w');
            for ct = 1:200
                if ct>1,tf.writeDirectory;end
                tf.setTag('Photometric',Tiff.Photometric.MinIsBlack);
                tf.setTag('Compression',Tiff.Compression.None);
                tf.setTag('BitsPerSample',16);
                tf.setTag('SamplesPerPixel',1);
                tf.setTag('SampleFormat',Tiff.SampleFormat.UInt);
                tf.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);
                tf.setTag('ImageLength',size(IM,1));
                tf.setTag('ImageWidth',size(IM,2));
                tf.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
                tf.setTag('ImageDescription',sprintf('ImageJ=1.51j\nchannels=%.0f',size(IM,3)));
                tf.write(IM(:,:,ct));
                t(ct)=toc;
            end
            tf.close();
        case 'fTIF'
            %timing vector
            t = zeros(1,size(IM,3)+1);
            tic
            fTIF = Fast_Tiff(filename);
            for ct = 1:size(IM,3)
                fTIF = fTIF.WriteIMG(IM(:,:,ct)');
                t(ct)=toc;
            end
            tic
            fTIF.close;
            toc
        otherwise
            error('unknown method')
    end
    S = (size(IM,1)*size(IM,2)*2)/2^20; %MB/frame
    y = S./diff(t);
    subplot(1,length(methods),M)
    plot([1:length(y)],y);
    title(sprintf('Writing with %s; mean = %.2f MB/s',method,mean(y)))
    ylabel('Writing speed (MB/s)')
    xlabel('Frame');
    drawnow;
end
