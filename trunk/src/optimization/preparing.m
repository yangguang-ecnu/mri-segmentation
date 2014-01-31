function [in out] = preparing

global source_tri
global var_array1
global sag_m1
global axial_m1
global vol_ax
global vol_sag
global vol_cor
global sub_1
global sub_2
global sub_3
global t
global current_tr
global c2b_coord
global list_edges
global neig_ax
global neig_sag
global neig_cor
global gradx_ax
global grady_ax
global gradx_sag
global grady_sag
global gradx_cor
global grady_cor

% Define the input struct
in.size_source = size(source_tri.X,1);

in.w1 = zeros(4,3*size(source_tri.X,1));
in.w2 = zeros(4,3*size(source_tri.X,1));

% in.w1(4,:) = 1;
% in.w2(4,:) = 1;

in.rows = size(vol_ax,1);
in.cols = size(vol_ax,2);

in.J1 = zeros(size(var_array1,1),3*size(source_tri.X,1));
in.J2 = zeros(size(var_array1,1),3*size(source_tri.X,1));

in.n = size(var_array1,1);

in.source_tri = source_tri;
in.var_array1 = var_array1;
in.sag_m1 = sag_m1;
in.axial_m1 = axial_m1;

for i = 1:in.n
    in.axial_m1_par{i} = in.axial_m1{sub_3(i)};
    in.sag_m1_par{i}   = in.sag_m1{sub_2(i)};
end

in.vol_ax = vol_ax;
in.vol_sag = vol_sag;

in.gradx_ax  = gradx_ax;
in.grady_ax  = grady_ax;
in.gradx_sag = gradx_sag;
in.grady_sag = grady_sag;
in.gradx_cor = gradx_cor;
in.grady_cor = grady_cor;

in.sub_1 = sub_1;
in.sub_2 = sub_2;
in.sub_3 = sub_3;
in.t = t;
in.current_tr = current_tr;
in.c2b_coord = c2b_coord;
in.list_edges = list_edges;

in.lambda = .01;
in.sqrt_lambda = sqrt(in.lambda);

for i = 1:length(list_edges)
    in.scalar(i) = -1/length(list_edges{i});
end

%% Neighborhood
size_neig = [5 5];

for k = 1:size(vol_ax,3);
    neig_ax{k}  = colfilt( vol_ax(:,:,k), size_neig,'sliding',@mean);
end

for k = 1:size(vol_sag,3);
    neig_sag{k} = colfilt( vol_sag(:,:,k),size_neig,'sliding',@mean);
end

for k = 1:size(vol_cor,3);
    neig_cor{k} = colfilt( vol_cor(:,:,k),size_neig,'sliding',@mean);
end

in.neig_ax  = neig_ax;
in.neig_sag = neig_sag;
in.neig_cor = neig_cor;


% Define the output struct
out.F = zeros(size(var_array1,1) + size(source_tri.X,1) + size(source_tri.X,1),1);
out.J = zeros(size(var_array1,1) + size(source_tri.X,1) + size(source_tri.X,1),6*size(source_tri.X,1)); % the control points are duplicated [X0;X0]
% out.F = zeros(size(var_array1,1),1);
% out.J = zeros(size(var_array1,1),6*size(source_tri.X,1)); % the control points are duplicated [X0;X0]

