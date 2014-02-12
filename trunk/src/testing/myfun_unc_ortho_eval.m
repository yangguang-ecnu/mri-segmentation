function F_obj = myfun_unc_ortho_eval(mesh0_X)


global source_tri
global var_array1
global var_array2
global var_array3

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global sub_2_int1
global sub_3_int1

global sub_2_int2
global sub_3_int2

global sub_2_int3
global sub_3_int3
global current_tr_int1
global c2b_coord_int1
global current_tr_int2
global c2b_coord_int2
global current_tr_int3
global c2b_coord_int3
global list_edges

global axial_M1
global sag_M1
global cor_M1


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

F1  = zeros(size(var_array1,1),1); % + size(source_tri.X,1) + size(source_tri.X,1),1);
F2  = zeros(size(var_array2,1),1); 
F3  = zeros(size(var_array3,1),1); 
F_s = zeros(size(mesh0_X,1),1);

r_ax = size(vol_ax_eval,1);
c_ax = size(vol_ax_eval,2);

r_sag = size(vol_sag_eval,1);
c_sag = size(vol_sag_eval,2);

r_cor = size(vol_cor_eval,1);
c_cor = size(vol_cor_eval,2);

%% First intersection %%
for i = 1: n1
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if var_array1(i,1) ~= -Inf
        tmp_v1 = axial_M1{sub_3_int1(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_M1_par1{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),r_ax);
        min_max_r2 = min(max(cl2,1),r_ax);
        min_max_c  = min(max(fl,1),c_ax);
        min_max_c2 = min(max(cl,1),c_ax);
        
        neig = [vol_ax_eval(min_max_r, min_max_c, sub_3_int1(i))   vol_ax_eval(min_max_r, min_max_c2, sub_3_int1(i));...
            vol_ax_eval(min_max_r2,min_max_c, sub_3_int1(i))   vol_ax_eval(min_max_r2,min_max_c2, sub_3_int1(i))];
        
        new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = sag_M1{sub_2_int1(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_M1_par1{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),r_sag);
        min_max_r2 = min(max(cl2,1),r_sag);
        min_max_c  = min(max(fl,1), c_sag);
        min_max_c2 = min(max(cl,1), c_sag);
        
        neig = [vol_sag_eval(min_max_r, min_max_c,sub_2_int1(i))   vol_sag_eval(min_max_r, min_max_c2,sub_2_int1(i));...
                vol_sag_eval(min_max_r2,min_max_c,sub_2_int1(i))   vol_sag_eval(min_max_r2,min_max_c2,sub_2_int1(i))];
        
        new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        %% Function
        F1(i) = new_im_ax - new_im_sag;
    else
        'hola1'
        F1(i) = 0;
    end
    
end

%% Second intersection %%
for i = 1: n2
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if var_array2(i,1) ~= -Inf
        tmp_v1 = axial_M1{sub_3_int2(i)} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_M1_par2{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),r_ax);
        min_max_r2 = min(max(cl2,1),r_ax);
        min_max_c  = min(max(fl,1),c_ax);
        min_max_c2 = min(max(cl,1),c_ax);
        
        neig = [vol_ax_eval(min_max_r, min_max_c, sub_3_int2(i))   vol_ax_eval(min_max_r, min_max_c2, sub_3_int2(i));...
                vol_ax_eval(min_max_r2,min_max_c, sub_3_int2(i))   vol_ax_eval(min_max_r2,min_max_c2, sub_3_int2(i))];
        
        new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 =  cor_M1{sub_2_int2(i)} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),r_cor);
        min_max_r2 = min(max(cl2,1),r_cor);
        min_max_c  = min(max(fl,1), c_cor);
        min_max_c2 = min(max(cl,1), c_cor);
        
        neig = [vol_cor_eval(min_max_r, min_max_c,sub_2_int2(i))   vol_cor_eval(min_max_r, min_max_c2,sub_2_int2(i));...
                vol_cor_eval(min_max_r2,min_max_c,sub_2_int2(i))   vol_cor_eval(min_max_r2,min_max_c2,sub_2_int2(i))];
        
        new_im_cor = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        %% Function
        F2(i) = new_im_ax - new_im_cor;
    else
        'hola2'
        F2(i) = 0;
    end
    
end
%% Third intersection %%
for i = 1: n3
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if var_array3(i,1) ~= -Inf
        tmp_v1 = sag_M1{sub_2_int3(i)} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_M1_par2{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),r_sag);
        min_max_r2 = min(max(cl2,1),r_sag);
        min_max_c  = min(max(fl,1),c_sag);
        min_max_c2 = min(max(cl,1),c_sag);
        
        neig = [vol_sag_eval(min_max_r, min_max_c, sub_2_int3(i))   vol_sag_eval(min_max_r, min_max_c2, sub_2_int3(i));...
            vol_sag_eval(min_max_r2,min_max_c, sub_2_int3(i))   vol_sag_eval(min_max_r2,min_max_c2, sub_2_int3(i))];
        
        new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = cor_M1{sub_3_int3(i)} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_M1_par2{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),r_cor);
        min_max_r2 = min(max(cl2,1),r_cor);
        min_max_c  = min(max(fl,1), c_cor);
        min_max_c2 = min(max(cl,1), c_cor);
        
        neig = [vol_cor_eval(min_max_r, min_max_c,sub_3_int3(i))   vol_cor_eval(min_max_r, min_max_c2,sub_3_int3(i));...
            vol_cor_eval(min_max_r2,min_max_c,sub_3_int3(i))   vol_cor_eval(min_max_r2,min_max_c2,sub_3_int3(i))];
        
        new_im_cor = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        %% Function
        F3(i) = new_im_sag - new_im_cor;
    else
        'hola3'
        F3(i) = 0; 
    end
    
end

%% Smooth term %%
for i = 1:size( source_tri.X,1)
    
    F_s( i )         =  sqrt_lambda * (sum( (tri1.X(i,:) - mean(tri1.X( list_edges{i},:))).^2 ));
    F_s(i + incr)    =  sqrt_lambda * (sum( (tri2.X(i,:) - mean(tri2.X( list_edges{i},:))).^2 ));
    F_s(i + 2*incr)  =  sqrt_lambda * (sum( (tri3.X(i,:) - mean(tri3.X( list_edges{i},:))).^2 ));
    
end


F_obj = sum(F1.^2) + sum(F2.^2) + sum(F3.^2) + sum(F_s.^2);

