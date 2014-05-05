clc;
close all;

%% Change the name of the file for saving the results
save_name = 'lassalas.mat';

%% Set up the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Set up parameters -------')
disp('-----------------------------------')

global optimizer

%% Set up the number of points per intersectin 't', 
%% 'lambda' for the regularization term
%% 'nx, ny, nz' for the grid dimensions
optimizer.t     = 54;
optimizer.lambda = .01;
optimizer.nxyz  = [5 5 5];

if ~exist('set_up','var')
    set_parameters;
end

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Start Optimization  -----')
disp('-----------------------------------')

%% Using gradient descent
optimizer.steepestdescent = 1;
optimizer.grad       = 1; %% 0/1 without gradient/provide gradient
optimizer.maxiter    = 40; %% maximum of iterations for the optimzer
optimizer.checkderiv = 'off'; %% 'on'/'off' check or not the provided gradient with MATLAB one
optimizer.tolfun     = 0.1; 
optimizer.tolx       = 0.01; 
optimizer.con        = 0; %% 0/1 -> fminunc / fmincon

%% Using metaheuristics, differential evolution
optimizer.metaheuristic = 0;

%% Set the number of views (default 3)
optimizer.views         = 3; %% provide the number of views used for the optimization, 2 or 3

start_optimization( );

%% Calculate the new images from the optimization result

disp('--------- Compute the new images and visualization --------------')
compute_new_images( );

%% Calculate the measurements
disp('--------- Compute error measurements --------------')
measurements(  );

%% Save the output if the name is provided
save(save_name, 'source_tri', 'target_tri_ax' , 'target_tri_sag', 'target_tri_cor',...
                'axial_m','axial_m1','sag_m','sag_m1','cor_m','cor_m1', ...
                'new_axial', 'new_sagittal', 'new_coronal', ...
                'vol_ax', 'vol_sag' , 'vol_cor',...
                'optimizer', 't_reg', 'measure_errors');
            
%% Choose the images to visualize, i.e, [9 8 8] corresponds to slice 9 of the axial view,
%% slice 8 from sagittal view and slice 8 from coronal one.
show_slices = [19 9 12];
visualization( views, show_slices );


