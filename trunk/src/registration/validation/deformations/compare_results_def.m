%% Main evaluation, compare the different points %%

clc
serie = 7;

lava_flex_n      = cellimages2mat(dcmdir.dcmPatient.Study.Series(serie,1).Images);
lava_flex_info = dcmdir.dcmPatient.Study.Series(serie,1).ImagesInfo;

[size1, size2, size3] = size(lava_flex_n);


%% deformation magnitude %%
sigma = 5;

% root_def = 'deformations/deformation_';
% root_res = 'deformations_results/result_deformation_';
root_def = 'Noise/noise_b';
root_res = 'Noise_results/deform_noiseb';

n_x = 4;
n_y = 4;
n_z = 3;

for i=1:length(sigma)
    
    disp(['--------- Deformation ', num2str(sigma(i)),' ----------------------------------']);
    %% registered triangulation
    file_name_def = strcat(root_def, num2str(sigma(i)), '.mat');
    file_name_res = strcat(root_res, num2str(sigma(i)), '.mat');
    
    def = load(file_name_def);
    res = load(file_name_res);
    
    save_name = strcat('error_deformations/error_deformation_', num2str(sigma(i)), '.mat');
    
    %% original and perturbed triangulation
    st  = coordinates_M_coordinates(lava_flex_n, lava_flex_info);
    tri = coordinates_deformation_def(st, res.source_tri_v, def );
    
    %% 6 comparisons
    %% First ax to sag
    current_tr = tsearchn(def.deform_tri_ax.X,def.deform_tri_ax.Triangulation,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);
    b          = cartToBary(res.source_tri_v, current_tr,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(res.source_tri_v, current_tr, b); 
    
    b2          = cartToBary(def.deform_tri_sag, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(def.deform_tri_sag, current_tr, b2);     
    
    
    current_tr1 = tsearchn(res.source_tri_v.X,res.source_tri_v.Triangulation,[def.var_array1_v(:,1) def.var_array1_v(:,2) def.var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
    b1 = cartToBary(res.target_tri_ax,  current_tr1,[def.var_array1_v(:,1) def.var_array1_v(:,2) def.var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri
    xc1_ax2  = baryToCart(def.deform_tri_ax,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
    xc1_sag2 = baryToCart(def.deform_tri_sag, current_tr1, b1);
    
    b          = cartToBary(res.target_tri_ax, current_tr,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp3   = baryToCart(res.target_tri_ax, current_tr, b); 
    b2          = cartToBary(def.deform_tri_sag, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(def.deform_tri_sag, current_tr, b); % b2
    
    
    error1 = sum((tri.xc1_sag - xc1_sag2).^2); %tri.xc1_sag
    
    %% Second sag to ax
    current_tr = tsearchn(def.deform_tri_sag.X,def.deform_tri_sag.Triangulation,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);
    b          = cartToBary(res.source_tri_v, current_tr,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(res.source_tri_v, current_tr, b); 
    
    b2           = cartToBary(def.deform_tri_ax, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(def.deform_tri_ax, current_tr, b2);     
    
    current_tr1 = tsearchn(res.source_tri_v.X,res.source_tri_v.Triangulation,[def.var_array1_v(:,1) def.var_array1_v(:,2) def.var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
    b1 = cartToBary(res.target_tri_sag,  current_tr1,[def.var_array1_v(:,1) def.var_array1_v(:,2) def.var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri
    xc1_ax2  = baryToCart(def.deform_tri_ax,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
%     xc1_sag2 = baryToCart(def.deform_tri_sag, current_tr1, b1);
    
    
    b          = cartToBary(res.target_tri_sag, current_tr,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp3   = baryToCart(res.target_tri_sag, current_tr, b); 
    b2           = cartToBary(def.deform_tri_ax, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(def.deform_tri_ax, current_tr, b2);
    
    error2 = sum((tri.xc1_ax - xc1_ax2).^2);%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Third ax to cor
    current_tr = tsearchn(def.deform_tri_ax.X,def.deform_tri_ax.Triangulation,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);
    b          = cartToBary(res.source_tri_v, current_tr,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(res.source_tri_v, current_tr, b); 
    
    b2           = cartToBary(def.deform_tri_cor, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(def.deform_tri_cor, current_tr, b2);     
    
    current_tr1 = tsearchn(res.source_tri_v.X,res.source_tri_v.Triangulation,[def.var_array2_v(:,1) def.var_array2_v(:,2) def.var_array2_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
    b1 = cartToBary(res.target_tri_ax,  current_tr1,[def.var_array2_v(:,1) def.var_array2_v(:,2) def.var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri
    xc2_cor2  = baryToCart(def.deform_tri_cor,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
%     xc1_sag2 = baryToCart(def.deform_tri_sag, current_tr1, b1);
    
    b          = cartToBary(res.target_tri_ax, current_tr,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp3   = baryToCart(res.target_tri_ax, current_tr, b); 
    b2          = cartToBary(def.deform_tri_cor, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(def.deform_tri_cor, current_tr, b2);
    
    
    error3 = sum((tri.xc2_cor - xc2_cor2).^2);%
    
    %% Fourth cor to ax
    current_tr = tsearchn(def.deform_tri_cor.X,def.deform_tri_cor.Triangulation,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);
    b          = cartToBary(res.source_tri_v, current_tr,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(res.source_tri_v, current_tr, b); 
    
    b2           = cartToBary(def.deform_tri_ax, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(def.deform_tri_ax, current_tr, b2);     
    
    current_tr1 = tsearchn(res.source_tri_v.X,res.source_tri_v.Triangulation,[def.var_array2_v(:,1) def.var_array2_v(:,2) def.var_array2_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
    b1 = cartToBary(res.target_tri_cor,  current_tr1,[def.var_array2_v(:,1) def.var_array2_v(:,2) def.var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri
    xc2_ax2  = baryToCart(def.deform_tri_ax,  current_tr1, b1); % cartesian codinates of the points wrt target tri ( p' )
%     xc1_sag2 = baryToCart(def.deform_tri_sag, current_tr1, b1);
    
    b          = cartToBary(res.target_tri_cor, current_tr,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp3   = baryToCart(res.target_tri_cor, current_tr, b); 
    b2           = cartToBary(def.deform_tri_ax, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(def.deform_tri_ax, current_tr, b2);
    
    error4 = sum((tri.xc2_ax - xc2_ax2).^2);% tri.xc2_ax
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Third sag to cor
    current_tr = tsearchn(def.deform_tri_sag.X,def.deform_tri_ax.Triangulation,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);
    b          = cartToBary(res.source_tri_v, current_tr,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(res.source_tri_v, current_tr, b); 
    
    b2          = cartToBary(def.deform_tri_cor, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(def.deform_tri_cor, current_tr, b2);     
    
    current_tr1 = tsearchn(res.source_tri_v.X,res.source_tri_v.Triangulation,[def.var_array3_v(:,1) def.var_array3_v(:,2) def.var_array3_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
    b1 = cartToBary(res.target_tri_sag,  current_tr1,[def.var_array3_v(:,1) def.var_array3_v(:,2) def.var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri
    xc3_cor2  = baryToCart(def.deform_tri_cor,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
%     xc1_sag2 = baryToCart(def.deform_tri_sag, current_tr1, b1);
    
    b          = cartToBary(res.target_tri_sag, current_tr,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp3   = baryToCart(res.target_tri_sag, current_tr, b); 
    b2          = cartToBary(def.deform_tri_cor, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(def.deform_tri_cor, current_tr, b2);
    
    
    error5 = sum((tri.xc3_cor - xc3_cor2).^2);% 
    
    %% Fourth cor to sag
    current_tr = tsearchn(def.deform_tri_cor.X,def.deform_tri_cor.Triangulation,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);
    b          = cartToBary(res.source_tri_v, current_tr,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(res.source_tri_v, current_tr, b); 
    
    b2           = cartToBary(def.deform_tri_sag, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(def.deform_tri_sag, current_tr, b2);     
    
    current_tr1 = tsearchn(res.source_tri_v.X,res.source_tri_v.Triangulation,[def.var_array3_v(:,1) def.var_array3_v(:,2) def.var_array3_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
    b1 = cartToBary(res.target_tri_cor,  current_tr1,[def.var_array3_v(:,1) def.var_array3_v(:,2) def.var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri
    xc3_sag2  = baryToCart(def.deform_tri_sag,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
%     xc1_sag2 = baryToCart(def.deform_tri_sag, current_tr1, b1);
    
    b           = cartToBary(res.target_tri_cor, current_tr,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp3   = baryToCart(res.target_tri_cor, current_tr, b); 
    b2          = cartToBary(def.deform_tri_sag, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(def.deform_tri_sag, current_tr, b2);
    
    error6 = sum((tri.xc3_sag - xc3_sag2).^2); %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    error = sum(error1 + error2 + error3 + error4 + error5 + error6)/(2 * size(tri.xc1_ax,1) + 2* size(tri.xc2_ax,1) + 2* size(tri.xc3_sag,1));
    save(save_name,'error');
    disp(['Error: ', num2str(error)]);
end
