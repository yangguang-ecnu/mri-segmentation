%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Imaging we have the deformation (target control points). 
%% Question: how to determine the location of a point p after the deformation,
%% p is not a source control point.
%% Answer: with barycentric coordinates
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;

global source_tri
global var_array
global axial_m1
global var_cell

rows = size(var_cell,1);
cols = size(var_cell,2);

tri = source_tri.Triangulation;

%% Plot the original triangulation and reference points.
figure
subplot(1,2,1);
tetramesh(source_tri); hold on;
alpha(.1)
plot3(var_array(:,1), var_array(:,2),var_array(:,3),'*r'); 
hold off;
axis equal;

%% Stretch the triangulation and compute the mapped locations of the incenters on the deformed triangulation.
% noise_data = source_tri.X + randn(size(source_tri.X));

target_tri = TriRep(tri,source_tri.X(:,1), source_tri.X(:,2),source_tri.X(:,3));
% target_tri = TriRep(tri,noise_data(:,1),noise_data(:,2),noise_data(:,3));

for i = 1:size(var_array,1)
    current_tr = tsearchn(source_tri.X,source_tri.Triangulation,var_array(i,:));
    %current_tr = tsearchn(noise_data,source_tri.Triangulation,var_array(i,:));
    v = 1:length(current_tr);
    b{i} = cartToBary(source_tri,v',var_array(i,:)); % barycentric coordinates of the points wrt source tri
    xc{i} = baryToCart(target_tri, v', b{i}); % cartesian coordinates of the points wrt target tri ( p' )
end

%% Get the intensity values of the p' in the first direction

for i = 1:size(var_array,1)
    
    [sub_tmp1 sub_tmp2 sub_tmp3] = ind2sub([size(var_cell,3) size(var_cell,2) size(var_cell,1)],i);
    
    tmp_v1 = axial_m1{sub_tmp1} * [xc{i} 1]'; % 3D point to 2D point in the frame coordinates

    tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)

    %% Make sure the indexes are correct and use bilinear interpolation

    neig = [vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(floor(tmp_v(2) + 1),1),cols),sub_tmp1) vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(ceil(tmp_v(2) + 1),1),cols),sub_tmp1);...
        vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(floor(tmp_v(2) + 1),1),cols),sub_tmp1) vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(ceil(tmp_v(2) + 1),1),cols),sub_tmp1)];
    bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
end            
% Plot the deformed triangulation and mapped locations of the reference points.
subplot(1,2,2);
tetramesh(target_tri); 
alpha(.1)
hold on;
for i = 1:size(var_array,1)
    plot3(xc{i}(:,1), xc{i}(:,2), xc{i}(:,3), '*r'); 
end
hold off;
axis equal;