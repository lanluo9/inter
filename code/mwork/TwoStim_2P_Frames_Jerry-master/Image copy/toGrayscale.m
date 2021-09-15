clear;
for i = 6 : 7 %loop through images that need to be changed
    filepath = strcat('coloredImages/', int2str(i), '.jpg'); %creates filepath for location of colored image. will need to create "Images" folder if doesn't exist already
    I = rgb2gray(imread(filepath));
    imshow(I)
    imsave;
end

%to make the grayscale images work with mworks, change the .j2c extension to .jpg