function [out] = coordinates_deformation_def(st, control_tri, def)

%% Calculate the intersection between the different planes

options = optimset('Display','off');

ax  = st.slices_ax;
sag = st.slices_sag;
cor = st.slices_cor;


current_tr1 = tsearchn(control_tri.X,control_tri.Triangulation,[def.var_array1_v(:,1) def.var_array1_v(:,2) def.var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
current_tr2 = tsearchn(control_tri.X,control_tri.Triangulation,[def.var_array2_v(:,1) def.var_array2_v(:,2) def.var_array2_v(:,3)]);
current_tr3 = tsearchn(control_tri.X,control_tri.Triangulation,[def.var_array3_v(:,1) def.var_array3_v(:,2) def.var_array3_v(:,3)]);

% New data positions with respect the new computed meshes
b1 = cartToBary(control_tri,  current_tr1,[def.var_array1_v(:,1) def.var_array1_v(:,2) def.var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc1_ax  = baryToCart(def.deform_tri_ax,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
xc1_sag = baryToCart(def.deform_tri_sag, current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )

b2 = cartToBary(control_tri,  current_tr2,[def.var_array2_v(:,1) def.var_array2_v(:,2) def.var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc2_ax  = baryToCart(def.deform_tri_ax, current_tr2, b2); % cartesian coordinates of the points wrt target tri ( p' )
xc2_cor = baryToCart(def.deform_tri_cor, current_tr2, b2); % cartesian coordinates of the points wrt target tri ( p' )

b3 = cartToBary(control_tri,  current_tr3,[def.var_array3_v(:,1) def.var_array3_v(:,2) def.var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc3_sag = baryToCart(def.deform_tri_sag, current_tr3, b3); % cartesian coordinates of the points wrt target tri ( p' )
xc3_cor = baryToCart(def.deform_tri_cor, current_tr3, b3); % cartesian coordinates of the points wrt target tri ( p' )

out.xc1_ax  = xc1_ax;
out.xc1_sag = xc1_sag;
out.xc2_ax  = xc2_ax;
out.xc2_cor = xc2_cor;
out.xc3_sag = xc3_sag;
out.xc3_cor = xc3_cor;



