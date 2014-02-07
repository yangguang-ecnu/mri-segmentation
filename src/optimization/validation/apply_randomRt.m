function out_plane = apply_randomRt(plane)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Apply a random rotation and translation to the given plane, define
%% by its 4 corners
%%
%% Inputs:  1. plane     -> 4 x 3 array, each row represents the 4 corners
%% Outputs: 1. out_plane -> 4 x 3 array, each row represents the 4 corners
%%                          of the new plane
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = -10;
b = 10;

[R, ~] = qr(randn(3));
t = a + (b-a).*rand(3,1);

for i = 1:size(plane,1)
    
    out_plane = R*plane(i,:)' + t;
    
end