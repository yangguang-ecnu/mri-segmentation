function out = add_rician_noise(im, perc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Add Rician noise to a given image, using the given sigma value
%% common noise in MRI data 
%%
%% Inputs: 1. im    -> m x n matrix (image)
%%         2. perc  -> percentage of Rician noise to be added 
%% Output: 1. out   -> m x n matrix (image) after applying noise
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = size(im);

% Create noisy data with Rician noise
level = perc * max(im(:))/100;
out = sqrt((im +level*randn(s)).^2+(level*randn(s)).^2);

