function out = add_rician_noise(im, perc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Add Rician noise to a given image, using the given sigma value
%% common noise in MRI data f(x;sigma) = 2 x-mu/sigma^2 x e^{-(x-mu)^2/(2sigma^2)}
%%
%% Inputs: 1. im    -> m x n matrix (image)
%%         2. mu    -> real value for  Rayleigh distribution 
%%         3. sigma -> real value for Rayleigh distribution  
%% Output: 1. out   -> m x n matrix (image) after applying noise
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r c] = size(im);
out = zeros(r,c);

perc = 2;
s=size(im);
% Create noisy data with Rician noise
level = perc * max(im(:))/100;
out = sqrt((im +level*randn(s)).^2+(level*randn(s)).^2);

