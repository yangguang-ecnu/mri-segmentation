function prepare_optimization( var_cell1_v, var_cell2_v, var_cell3_v, var_array1_v, var_array2_v, var_array3_v, source_tri_v )

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------- Start Optimization  -----')

disp('--------- Preparing first intersection: axial/sagittal  -----')
size1_int1 = size(var_cell1_v,3);
size2_int1 = size(var_cell1_v,2);
size3_int1 = size(var_cell1_v,1);

global sub_1_int1_v
global sub_2_int1_v
global sub_3_int1_v
global current_tr_int1_v
global c2b_coord_int1_v


[sub_1_int1_v, sub_2_int1_v, sub_3_int1_v] = ind2sub([size1_int1 size2_int1 size3_int1],1:size(var_array1_v,1));
current_tr_int1_v = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int1_v  = cartToBary(source_tri_v,current_tr_int1_v,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri

disp('--------- Preparing second intersection: axial/coronal  -----')
size1_int2 = size(var_cell2_v,3);
size2_int2 = size(var_cell2_v,2);
size3_int2 = size(var_cell2_v,1);

global sub_1_int2_v
global sub_2_int2_v
global sub_3_int2_v
global current_tr_int2_v
global c2b_coord_int2_v


[sub_1_int2_v, sub_2_int2_v, sub_3_int2_v] = ind2sub([size1_int2 size2_int2 size3_int2],1:size(var_array2_v,1));
current_tr_int2_v = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int2_v  = cartToBary(source_tri_v,current_tr_int2_v,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri


disp('--------- Preparing third intersection: sagittal/coronal  -----')
size1_int3 = size(var_cell3_v,3);
size2_int3 = size(var_cell3_v,2);
size3_int3 = size(var_cell3_v,1);

global sub_1_int3_v
global sub_2_int3_v
global sub_3_int3_v
global current_tr_int3_v
global c2b_coord_int3_v


[sub_1_int3_v, sub_2_int3_v, sub_3_int3_v] = ind2sub([size1_int3 size2_int3 size3_int3],1:size(var_array3_v,1));
current_tr_int3_v = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int3_v  = cartToBary(source_tri_v,current_tr_int3_v,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri


preparing_eval;

