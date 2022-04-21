% attempt to denoise tif for cellpose: did not work
net = denoisingNetwork('dncnn');
img = imread('cellpose_stim_resp_gauss.tif');
denoised_img = denoiseImage(img, net);

figure
imshowpair(img(1:100, 1:100)',denoised_img(1:100, 1:100)','montage');
set(gcf, 'Position', get(0, 'Screensize'));
