%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main reading patient
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = 'resources/patients/DALLOT/';

%% Read the dicomdir and save the info in dicomdir.mat in 'folder'
load_dicomdir(folder);

mat_file = strcat(folder,'dicomdir.mat');
load(mat_file);

%% Define a struct with axial, sagittal and coronal views
views = get_3views(dcmdir);

%% Clear the rest of the data

clear dcmdir dcmdirlistVal folder mat_file
clc