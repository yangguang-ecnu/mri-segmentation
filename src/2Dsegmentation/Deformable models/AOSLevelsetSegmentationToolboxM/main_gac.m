%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main program for GAC 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
%% Select an image
I = double(views.axial(:,:,17));
I = anisodiff2D(I, 10, 1/7, 50, 2);
%I = patient_or.sagittal(:,:,10);
%I = patient_or.coronal(:,:,12);
dims = size(I);

%% Set the initial level set
center = [249,189]; radius =  10;     
phi1 = ac_SDF_2D('circle', dims, center, radius);

%% Compute the gradient map 
gr = ac_gradient_map(I,.07);
figure;
imshow(gr,[]);

%% Set the GAC parameters
countour_weight = 5; expansion_weight = 1; 
delta_t = 5; n_iters = 100; show_result = 1; 

%% Apply GAC
phi = ac_GAC_model(gr, phi1, countour_weight, expansion_weight, ...
        delta_t, n_iters, show_result);
 
