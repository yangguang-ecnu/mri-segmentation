function [im,file]=dicom(filename,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Given a dicome file, it returns the image and its dicom info
%%
%% Inputs:  1. filename -> string with the name of the file
%%          2. show -> 1 if show image, 0 otherwise
%% Outputs: 1. im -> image
%%          2. file -> struct containing the dicom info
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file=dicominfo(filename);
im=dicomread(filename);

imval.Width = file.Width;
imval.Height = file.Height;
imval.BitDepth = file.BitDepth;
imval.PixelSpacing = file.PixelSpacing;
imval.BitsAllocated = file.BitsAllocated;
imval.BitsStored = file.BitsStored;
imval.HighBit = file.HighBit;
imval.PatientName = file.PatientName;

if nargin<2
    show = false;
end

if show
    figure;
    imshow(im,[]);
end