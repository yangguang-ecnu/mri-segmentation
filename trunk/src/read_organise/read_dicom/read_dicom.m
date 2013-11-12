function [images dicom_inf dcm] = read_dicom(folder,resize_factor,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read and show the DICOM images from
%% a given folder.
%% Inputs: 1. folder -> name of folder containing dcm files
%%         2. resize_factor -> resize images from dcm files
%%         3. show ->  0/1 = show/not show the images
%%
%% Output: 1. images -> images obtain from dcm files and resize
%%                      given the resize_factor
%%         2. dicom_inf -> dicom information
%%         3. dcm -> all files in folder
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dcm = dir(folder);

for i=3:size(dcm,1)
    
    [im inf] = dicom(dcm(i).name,0);
    im = imresize(im,resize_factor,'bicubic');
    images{i-2} = im;
    dicom_inf{i-2} = inf;

    %% Check if it is the same patien
    if show == 1
        imshow(im,[]);title(['Image: ',dcm(i).name]);
        pause;
    end
end


if show == 1
    close all;
end