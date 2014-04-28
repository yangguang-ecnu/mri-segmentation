clc;
close all;

save_name = 'giraud.mat';

%% Set up the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Set up parameters -------')
disp('-----------------------------------')

if ~exist('set_up','var')
    set_parameters;
end

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Start Optimization  -----')
disp('-----------------------------------')

global optimizer

optimizer.steepestdescent = 1;
optimizer.views      = 3; %% provide the number of views used for the optimization, 2 or 3
optimizer.grad       = 1; %% 0/1 without gradient/provide gradient
optimizer.maxiter    = 40; %% maximum of iterations for the optimzer
optimizer.checkderiv = 'off'; %% 'on'/'off' check or not the provided gradient with MATLAB one
optimizer.tolfun     = 0.01; 
optimizer.tolx       = 0.001; 
optimizer.con        = 0; %% 0/1 -> fminunc / fmincon

optimizer.metaheuristic = 0;
optimizer.views         = 3; %% provide the number of views used for the optimization, 2 or 3

start_optimization( );

%% Calculate the new images from the optimization result

disp('--------- Compute the new images and visualization --------------')
compute_new_images( );

%% Save the output if the name is provided

global target_tri_ax
global target_tri_sag
global target_tri_cor

global source_tri

global new_axial
global new_sagittal
global new_coronal

save(save_name, 'source_tri', 'target_tri_ax' , 'target_tri_sag', 'target_tri_cor',...
                'axial_m','axial_m1','sag_m','sag_m1','cor_m','cor_m1', ...
                'new_axial', 'new_sagittal', 'new_coronal', ...
                'vol_ax', 'vol_sag' , 'vol_cor');


%% Choose the images to visualize, i.e, [9 8 8] corresponds to slice 9 of the axial view,
%% slice 8 from sagittal view and slice 8 from coronal one.
visualization( views, [8 9 9] );





