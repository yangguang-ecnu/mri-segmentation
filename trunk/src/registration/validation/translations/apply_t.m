function [out_plane, tr] = apply_t(plane,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Apply a random rotation and translation to the given plane, define
%% by its 4 corners
%%
%% Inputs:  1. plane     -> 4 x 3 array, each row represents the 4 corners
%%          2. t         -> 3 X 1 array, the translation vector  
%% Outputs: 1. out_plane -> 4 x 3 array, each row represents the 4 corners
%%                          of the new plane
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(plane,1)
    
    out_plane(i,:) = plane(i,:)' + t;
    
end

tr = eye(4);
tr(1:3,4) = t;