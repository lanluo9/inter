
close all

exptDateYYMMDD = 230330;
objective = '16x';
zoom = 2.0;

cellpose_mask_mat = 'Z:\All_Staff\home\lan\Analysis\2P\230330_i1380\cellpose_mask.mat';
tmp = load(cellpose_mask_mat);
fov_img = tmp.cellpose_mask;

[x_um, y_um, sb_img_50um] = find_2p_fov_size_in_mm(exptDateYYMMDD, objective, fov_img, zoom);

%%

scalebar_50um_pixel = 77.2816 % printed from find_2p_fov_size_in_mm
fov_img_xaxis_pixel = size(fov_img, 2); % fov_img = [y pixel, x pixel]

fov_img_xaxis_inch = 11.639; % measured from corel draw, copied from cellpose mask i1380
scalebar_50um_inch = scalebar_50um_pixel / fov_img_xaxis_pixel * fov_img_xaxis_inch;

scalebar_10um_inch = scalebar_50um_inch / 5


%%