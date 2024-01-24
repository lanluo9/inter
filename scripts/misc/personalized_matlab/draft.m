% Observed data

%% vis p<0.05 bonf
% % V1 vs LM
% n1 = 286; N1 = 583;
% n2 = 620; N2 = 1715;

% % V1 vs LI
% n1 = 286; N1 = 583;
% n2 = 123; N2 = 691;

% LM vs LI
% n1 = 620; N1 = 1715;
% n2 = 123; N2 = 691;

%% well tuned isi=6k

% % % V1 vs LM
% n1 = 296; N1 = 449;
% n2 = 578; N2 = 983;

% % % V1 vs LI
% n1 = 296; N1 = 449;
% n2 = 96; N2 = 222;

% % LM vs LI
% n1 = 578; N1 = 983;
% n2 = 96; N2 = 222;

% %% well tuned isi=250
% % V1 vs LM
% n1 = 281; N1 = 449;
% n2 = 532; N2 = 983;

% % % V1 vs LI
% n1 = 281; N1 = 449;
% n2 = 112; N2 = 222;

% % LM vs LI
% n1 = 532; N1 = 983;
% n2 = 112; N2 = 222;

x1 = [repmat('a',N1,1); repmat('b',N2,1)];

x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];

[tbl, chi2stat, pval] = crosstab(x1,x2)



% 
% close all
% 
% exptDateYYMMDD = 230330;
% objective = '16x';
% zoom = 2.0;
% 
% cellpose_mask_mat = 'Z:\All_Staff\home\lan\Analysis\2P\230330_i1380\cellpose_mask.mat';
% tmp = load(cellpose_mask_mat);
% fov_img = tmp.cellpose_mask;
% 
% [x_um, y_um, sb_img_50um] = find_2p_fov_size_in_mm(exptDateYYMMDD, objective, fov_img, zoom);
% 
% %%
% 
% scalebar_50um_pixel = 77.2816 % printed from find_2p_fov_size_in_mm
% fov_img_xaxis_pixel = size(fov_img, 2); % fov_img = [y pixel, x pixel]
% 
% fov_img_xaxis_inch = 11.639; % measured from corel draw, copied from cellpose mask i1380
% scalebar_50um_inch = scalebar_50um_pixel / fov_img_xaxis_pixel * fov_img_xaxis_inch;
% 
% scalebar_10um_inch = scalebar_50um_inch / 5
% 
% 
% 