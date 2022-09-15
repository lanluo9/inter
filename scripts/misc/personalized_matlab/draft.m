% % attempt to denoise tif for cellpose: did not work
% net = denoisingNetwork('dncnn');
% img = imread('cellpose_stim_resp_gauss.tif');
% denoised_img = denoiseImage(img, net);
% 
% figure
% imshowpair(img(1:100, 1:100)',denoised_img(1:100, 1:100)','montage');
% set(gcf, 'Position', get(0, 'Screensize'));


% test large sbx
% t = rand(264, 796, 100000, 3); % cannot leave 3 or 4 of sbx sized matrix
% in memory

t = rand(170, 100000, 3);
