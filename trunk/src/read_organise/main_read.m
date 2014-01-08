%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main reading patient
%%
%% Execute:
%% -change the folder path (line 10) where the patient is. 
%%  This folder should contain the dicomdir file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
folder = 'resources/patients/SAPHO/';%'/media/IMATION HDD/PhD_CF/backupAxis1/MRIs/ICE 2 - 2011/fibrome/';

%% Read the dicomdir and save the info in dicomdir.mat in 'folder'
load_dicomdir(folder);

mat_file = strcat(folder,'dicomdir.mat');
load(mat_file);

%% Define a struct with axial, sagittal and coronal views
views = get_3views(dcmdir);

%% Clear the rest of the data

clear dcmdirlistVal folder mat_file %dcmdir
clc

