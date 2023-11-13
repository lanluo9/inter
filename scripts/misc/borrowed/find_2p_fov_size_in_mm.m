function [x_um, y_um, sb_img_50um] = find_2p_fov_size_in_mm(exptDateYYMMDD, obj, img, zm)%% for 30 Hz imaging on 2P

% % originally scalebarCalib.m
% % commented and modified by LL, 2023-11

% % usage: this function makes a scale bar for two photon field of view (2p FOV)

% % input:
% exptDateYYMMDD (experiment date, YYMMDD)
% obj (objective, like '16x'), 
% img (field of view as an image, of size [x_pixel, y_pixel])
% zm (zoom, like 1.7, 2, etc)

% % output:
% x_um, y_um (unsure what this is)
% sb_img_50um (probably scale bar of image, equivalent to real life 50 um)

d = datestr(datenum(num2str(exptDateYYMMDD,'%d'),'yymmdd'),'yymmdd');
dprime = str2num(d(1:2));
if strcmp(obj,'16x')
    if dprime<=16
        % 16x 2014 - 2016
        x_um = 555;
        y_um = 233;
    elseif dprime>16
        % 16x 2017-2019
        x_um = 1030;
        y_um = 581;
    end
elseif strcmp(obj,'25x')
    if dprime<=15
        error('scale info not available for this date')
    elseif dprime>16
        % 25x 2017-2019
        x_um = 673; 
        y_um = 382;
    end
end
if length(nargin > 3)
    x_um = x_um/zm;
    y_um = y_um/zm;
end

%%
if nargin > 2
    [y_pix,x_pix] = size(img);
    x_calib = x_pix/x_um;
    y_calib = y_pix/y_um;
    x50um = x_calib*50;
    y5um = y_calib*5;

    sb_x_ind = 500:500+x50um;
    sb_y_ind = 200:200+y5um;

    sb_img = zeros(y_pix,x_pix);
    sb_img(sb_y_ind,sb_x_ind) = 1;
    sb_img_50um = sb_img;
%     figure;imagesc(sb_img)
end