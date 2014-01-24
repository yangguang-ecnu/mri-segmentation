function [in out] = preparing

global source_tri
global var_array1
global sag_m1
global axial_m1
global vol_ax
global vol_sag
global sub_1
global sub_2
global sub_3
global t
global current_tr
global c2b_coord
global list_edges

% Define the input struct
in.size_source = size(source_tri.X,1);

in.w1 = zeros(3,3*size(source_tri.X,1));
in.w2 = zeros(3,3*size(source_tri.X,1));

in.rows = size(vol_ax,1);
in.cols = size(vol_ax,2);

in.J1 = zeros(size(var_array1,1),3*size(source_tri.X,1));
in.J2 = zeros(size(var_array1,1),3*size(source_tri.X,1));

in.n = size(var_array1,1);

in.source_tri = source_tri;
in.var_array1 = var_array1;
in.sag_m1 = sag_m1;
in.axial_m1 = axial_m1;
in.vol_ax = vol_ax;
in.vol_sag = vol_sag;
in.sub_1 = sub_1;
in.sub_2 = sub_2;
in.sub_3 = sub_3;
in.t = t;
in.current_tr = current_tr;
in.c2b_coord = c2b_coord;
in.list_edges = list_edges;

for i = 1:length(list_edges)
    in.scalar(i) = -1/length(list_edges{i});
end
in.lambda = 1;


% Define the output struct
out.F = zeros(size(var_array1,1) + size(source_tri.X,1),1);
out.J = zeros(size(var_array1,1) + size(source_tri.X,1),6*size(source_tri.X,1)); % the control points are duplicated [X0;X0]

