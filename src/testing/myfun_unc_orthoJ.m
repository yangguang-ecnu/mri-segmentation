function [F_obj J] = myfun_unc_orthoJ(mesh0_X)


global source_tri
global var_array1
global var_array2
global var_array3
% global sag_m1
% global axial_m1
% global cor_m1
global vol_ax
global vol_sag
global vol_cor
global sub_1_int1
global sub_2_int1
global sub_3_int1
global sub_1_int2
global sub_2_int2
global sub_3_int2
global sub_1_int3
global sub_2_int3
global sub_3_int3
global current_tr_int1
global c2b_coord_int1
global current_tr_int2
global c2b_coord_int2
global current_tr_int3
global c2b_coord_int3
global list_edges
global gradx_ax
global gradx_sag
global gradx_cor
global grady_ax
global grady_sag
global grady_cor
global axial_m1
global sag_m1
global cor_m1
global scalar
% global axial_m1_par1
% global sag_m1_par1
% global cor_m1_par1
% 
% global axial_m1_par2
% global sag_m1_par2
% global cor_m1_par2

tri1 = TriRep( source_tri.Triangulation,[mesh0_X(1:size( source_tri.X,1),:) source_tri.X(:,3)]); % define the new mesh for the axial
tri2 = TriRep( source_tri.Triangulation,[source_tri.X(:,1) mesh0_X(size( source_tri.X,1)+1:2 * size( source_tri.X,1),:)]); % define the new mesh for the sagittal
tri3 = TriRep( source_tri.Triangulation,[mesh0_X(2 * size( source_tri.X,1)+1:end,1) source_tri.X(:,2) mesh0_X(2 * size( source_tri.X,1)+1:end,2)]); % define the new mesh for the sagittal

% First intersection
b2c_ncoord1_int1 = baryToCart(tri1,  current_tr_int1,  c2b_coord_int1);
b2c_ncoord2_int1 = baryToCart(tri2,  current_tr_int1,  c2b_coord_int1);

% Second intersection
b2c_ncoord1_int2 = baryToCart(tri1,  current_tr_int2,  c2b_coord_int2);
b2c_ncoord2_int2 = baryToCart(tri3,  current_tr_int2,  c2b_coord_int2);

% Third intersection
b2c_ncoord1_int3 = baryToCart(tri2,  current_tr_int3,  c2b_coord_int3);
b2c_ncoord2_int3 = baryToCart(tri3,  current_tr_int3,  c2b_coord_int3);

lambda = .01;
sqrt_lambda = sqrt(lambda);

n1 = size(var_array1,1);
n2 = size(var_array2,1);
n3 = size(var_array3,1);

incr = size(source_tri.X,1);

% incr_row = in.n + in.size_source;

F1  = zeros(size(var_array1,1),1); % + size(source_tri.X,1) + size(source_tri.X,1),1);
F2  = zeros(size(var_array2,1),1); 
F3  = zeros(size(var_array3,1),1); 
F_s = zeros(size(mesh0_X,1),1);

rows = size(vol_ax,1);
cols = size(vol_ax,2);

%% First intersection %%
for i = 1: n1
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = axial_m1{sub_3_int1(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par1{i}
    
    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r  = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c  = min(max(fl,1),cols);
    min_max_c2 = min(max(cl,1),cols);
    
    neig = [vol_ax(min_max_r, min_max_c, sub_3_int1(i))   vol_ax(min_max_r, min_max_c2, sub_3_int1(i));...
            vol_ax(min_max_r2,min_max_c, sub_3_int1(i))   vol_ax(min_max_r2,min_max_c2, sub_3_int1(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    grad_ax  = [gradx_ax( min_max_r, min_max_c, sub_3_int1(i)) grady_ax( min_max_r, min_max_c, sub_3_int1(i))];
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = sag_m1{sub_2_int1(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1_par1{i}

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);

    min_max_r  = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c  = min(max(fl,1), cols);
    min_max_c2 = min(max(cl,1), cols);
    
    neig = [vol_sag(min_max_r, min_max_c,sub_2_int1(i))   vol_sag(min_max_r, min_max_c2,sub_2_int1(i));...
            vol_sag(min_max_r2,min_max_c,sub_2_int1(i))   vol_sag(min_max_r2,min_max_c2,sub_2_int1(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    grad_sag = [gradx_sag(min_max_r, min_max_c, sub_2_int1(i)) grady_sag(min_max_r, min_max_c, sub_2_int1(i))];
     
    %% Function
    F1(i) = new_im_ax - new_im_sag;
    
    %% Jacobian
    w  = zeros(3, 9 * incr);
    w2 = zeros(3, 9 * incr);
    
    ind_w = (source_tri.Triangulation(current_tr_int1(i),:)-1)*3 + 1;
    
    w( 1,ind_w )            = c2b_coord_int1(i,:);
    w2(1,ind_w + 3 * incr)  = c2b_coord_int1(i,:);
    
    w(2,:)  = circshift(w(1,:)',1);
    w(3,:)  = circshift(w(1,:)',2);
    w2(2,:) = circshift(w2(1,:)',1);
    w2(3,:) = circshift(w2(1,:)',2);
    
    Js_int1(i,:) = (grad_ax * axial_m1{sub_3_int1(i)}(1:2,1:3) * w) - (grad_sag * sag_m1{sub_2_int1(i)}(1:2,1:3) * w2);
end

%% Second intersection %%
for i = 1: n2
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = axial_m1{sub_3_int2(i)} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par2{i}
    
    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r  = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c  = min(max(fl,1),cols);
    min_max_c2 = min(max(cl,1),cols);
    
    neig = [vol_ax(min_max_r, min_max_c, sub_3_int2(i))   vol_ax(min_max_r, min_max_c2, sub_3_int2(i));...
            vol_ax(min_max_r2,min_max_c, sub_3_int2(i))   vol_ax(min_max_r2,min_max_c2, sub_3_int2(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    grad_ax  = [gradx_ax( min_max_r, min_max_c, sub_3_int2(i)) grady_ax( min_max_r, min_max_c, sub_3_int2(i))];
    
    %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 =  cor_m1{sub_2_int2(i)} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);

    min_max_r  = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c  = min(max(fl,1), cols);
    min_max_c2 = min(max(cl,1), cols);
    
    neig = [vol_cor(min_max_r, min_max_c,sub_2_int2(i))   vol_cor(min_max_r, min_max_c2,sub_2_int2(i));...
            vol_cor(min_max_r2,min_max_c,sub_2_int2(i))   vol_cor(min_max_r2,min_max_c2,sub_2_int2(i))];
    
    new_im_cor = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    grad_cor  = [gradx_cor( min_max_r, min_max_c, sub_2_int2(i)) grady_cor( min_max_r, min_max_c, sub_2_int2(i))];
    
    %% Function
    F2(i) = new_im_ax - new_im_cor;
    
    %% Jacobian
    w  = zeros(3, 9 * incr);
    w2 = zeros(3, 9 * incr);
    
    ind_w = (source_tri.Triangulation(current_tr_int2(i),:)-1)*3 + 1;
    
    w( 1,ind_w )               = c2b_coord_int2(i,:);
    w2(1,ind_w + 6 * incr) = c2b_coord_int2(i,:);
    
    w(2,:) = circshift(w(1,:)',1);
    w(3,:) = circshift(w(1,:)',2);
    w2(2,:) = circshift(w2(1,:)',1);
    w2(3,:) = circshift(w2(1,:)',2);
    
    Js_int2(i,:) = (grad_ax * axial_m1{sub_3_int2(i)}(1:2,1:3) * w) - (grad_cor * cor_m1{sub_2_int2(i)}(1:2,1:3) * w2);
    
end
%% Third intersection %%
for i = 1: n3
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = sag_m1{sub_2_int3(i)} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_m1_par2{i}
    
    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r  = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c  = min(max(fl,1),cols);
    min_max_c2 = min(max(cl,1),cols);
    
    neig = [vol_sag(min_max_r, min_max_c, sub_2_int3(i))   vol_sag(min_max_r, min_max_c2, sub_2_int3(i));...
            vol_sag(min_max_r2,min_max_c, sub_2_int3(i))   vol_sag(min_max_r2,min_max_c2, sub_2_int3(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    grad_sag  = [gradx_sag( min_max_r, min_max_c, sub_2_int3(i)) grady_sag( min_max_r, min_max_c, sub_2_int3(i))];
    
    %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = cor_m1{sub_3_int3(i)} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_m1_par2{i}

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);

    min_max_r  = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c  = min(max(fl,1), cols);
    min_max_c2 = min(max(cl,1), cols);
    
    neig = [vol_cor(min_max_r, min_max_c,sub_3_int3(i))   vol_cor(min_max_r, min_max_c2,sub_3_int3(i));...
            vol_cor(min_max_r2,min_max_c,sub_3_int3(i))   vol_cor(min_max_r2,min_max_c2,sub_3_int3(i))];
    
    new_im_cor = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    grad_cor  = [gradx_cor( min_max_r, min_max_c, sub_3_int3(i)) grady_cor( min_max_r, min_max_c, sub_3_int3(i))];
    %% Function
    F3(i) = new_im_sag - new_im_cor;
    
    %% Jacobian
    w  = zeros(3, 9 * incr);
    w2 = zeros(3, 9 * incr);
    
    ind_w = (source_tri.Triangulation(current_tr_int3(i),:)-1)*3 + 1;
    
    w( 1,ind_w + 3 * incr) = c2b_coord_int3(i,:);
    w2(1,ind_w + 6 * incr) = c2b_coord_int3(i,:);
    
    w(2,:) = circshift(w(1,:)',1);
    w(3,:) = circshift(w(1,:)',2);
    w2(2,:) = circshift(w2(1,:)',1);
    w2(3,:) = circshift(w2(1,:)',2);
    
    Js_int3(i,:) = (grad_sag * sag_m1{sub_2_int3(i)}(1:2,1:3) * w) - (grad_cor * cor_m1{sub_3_int3(i)}(1:2,1:3) * w2);
    
end

incr_row = n1 + n2 + n3 + incr;
incr_col = 3 * incr;

%% Smooth term %%
for i = 1:size( source_tri.X,1)
    
    mean_tmp1 = mean(tri1.X(list_edges{i},:));
    mean_tmp2 = mean(tri2.X(list_edges{i},:));
    mean_tmp3 = mean(tri3.X(list_edges{i},:));
    lapl_tri1 = tri1.X(i,:) - mean_tmp1;
    lapl_tri2 = tri2.X(i,:) - mean_tmp2;
    lapl_tri3 = tri3.X(i,:) - mean_tmp3;
    
    F_s( i )         =  sqrt_lambda * (sum( lapl_tri1.^2 ));
    F_s(i + incr)    =  sqrt_lambda * (sum( lapl_tri2.^2 ));
    F_s(i + 2*incr)  =  sqrt_lambda * (sum( lapl_tri3.^2 ));
    

    jsi1 = (i - 1) * 3 + 1;
    jsi2 = jsi1 + incr_col;
    jsi3 = jsi2 + incr_col;
    
    Js_s(i    , jsi1:jsi1+2)        = (2 * sqrt_lambda).* lapl_tri1;%(in.sqrt_lambda * ( 1 / norm_lapl_tri1 ) ) .* lapl_tri1; 
    Js_s(i + incr, jsi2:jsi2+2)     = (2 * sqrt_lambda).* lapl_tri2;%(in.sqrt_lambda * ( 1 / norm_lapl_tri2 ) ) .* lapl_tri2; 
    Js_s(i + 2 * incr, jsi3:jsi3+2) = (2 * sqrt_lambda).* lapl_tri3;
    
    for j = 1:length(list_edges{i})
        edgein1 = (list_edges{i}(j)-1) * 3 + 1;
        edgein2 = edgein1 + incr_col;
        edgein3 = edgein2 + incr_col;
        Js_s(i   ,     edgein1:edgein1+2)     = Js_s(i ,     jsi1:jsi1+2).*[scalar(i) scalar(i) scalar(i)];
        Js_s(i + incr, edgein2:edgein2+2)     = Js_s(i + incr, jsi2:jsi2+2).*[scalar(i) scalar(i) scalar(i)];
        Js_s(i + 2 * incr, edgein3:edgein3+2) = Js_s(i + 2*incr, jsi3:jsi3+2).*[scalar(i) scalar(i) scalar(i)];
    end
    
end


F_obj = sum(F1.^2) + sum(F2.^2) + sum(F3.^2) + sum(F_s.^2);
J = sum(Js_int1) + sum(Js_int2) + sum(Js_int3) + sum(Js_s); %[Js_int1;Js_int2;Js_int3;Js_s];
