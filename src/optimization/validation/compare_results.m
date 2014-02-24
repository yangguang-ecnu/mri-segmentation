%% Main evaluation, compare the different points %%

clc
serie = 7;

lava_flex_n      = cellimages2mat(dcmdir.dcmPatient.Study.Series(serie,1).Images);
lava_flex_info = dcmdir.dcmPatient.Study.Series(serie,1).ImagesInfo;

[size1, size2, size3] = size(lava_flex_n);


%% deformation magnitude %%
sigma = [0 2 3 4 5 6 7 8];
% sigma = [1 2 3 4 5 6 7 8 9 10 11 15]; % noise
% 
t = 0:1/14:1;

root = 'deformation_';
% root = 'noise_';

n_x = 4;
n_y = 4;
n_z = 3;

for i=1%:length(sigma)
    
    disp(['--------- Deformation ', num2str(sigma(i)),' ----------------------------------']);
    %% registered triangulation
    file_name = strcat(root, num2str(sigma(i)), '.mat');
    reg = load(file_name);
    save_name = strcat('error_deform_', num2str(sigma(i)), '.mat');
    
    %% original and perturbed triangulation
    st  = coordinates_M_coordinates(lava_flex_n, lava_flex_info);
    tri = coordinates_deformation(t, st, n_x, n_y, n_z, sigma(i));
    
    %% 6 comparisons
    %% First ax to sag
    current_tr = tsearchn(tri.deform_tri_ax.X,tri.deform_tri_ax.Triangulation,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);
    b          = cartToBary(tri.control_tri, current_tr,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(tri.control_tri, current_tr, b); 
    
    b2           = cartToBary(tri.deform_tri_sag, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(tri.deform_tri_sag, current_tr, b2);     
    
    
    sag_comp3   = baryToCart(reg.target_tri_ax, current_tr, b); 
    b2           = cartToBary(tri.deform_tri_sag, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(tri.deform_tri_sag, current_tr, b2);
    
    
    error1 = sum((sag_comp2 - sag_comp4).^2);
    
    %% Second sag to ax
    current_tr = tsearchn(tri.deform_tri_sag.X,tri.deform_tri_sag.Triangulation,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);
    b          = cartToBary(tri.control_tri, current_tr,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(tri.control_tri, current_tr, b); 
    
    b2           = cartToBary(tri.deform_tri_ax, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(tri.deform_tri_ax, current_tr, b2);     
    
    
    sag_comp3   = baryToCart(reg.target_tri_sag, current_tr, b); 
    b2           = cartToBary(tri.deform_tri_ax, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(tri.deform_tri_ax, current_tr, b2);
    
    error2 = sum((sag_comp2 - sag_comp4).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Third ax to cor
    current_tr = tsearchn(tri.deform_tri_ax.X,tri.deform_tri_ax.Triangulation,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);
    b          = cartToBary(tri.control_tri, current_tr,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(tri.control_tri, current_tr, b); 
    
    b2           = cartToBary(tri.deform_tri_cor, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(tri.deform_tri_cor, current_tr, b2);     
    
    
    sag_comp3   = baryToCart(reg.target_tri_ax, current_tr, b); 
    b2          = cartToBary(tri.deform_tri_cor, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(tri.deform_tri_cor, current_tr, b2);
    
    
    error3 = sum((sag_comp2 - sag_comp4).^2);
    
    %% Fourth cor to ax
    current_tr = tsearchn(tri.deform_tri_cor.X,tri.deform_tri_cor.Triangulation,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);
    b          = cartToBary(tri.control_tri, current_tr,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(tri.control_tri, current_tr, b); 
    
    b2           = cartToBary(tri.deform_tri_ax, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(tri.deform_tri_ax, current_tr, b2);     
    
    
    sag_comp3   = baryToCart(reg.target_tri_cor, current_tr, b); 
    b2           = cartToBary(tri.deform_tri_ax, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(tri.deform_tri_ax, current_tr, b2);
    
    error4 = sum((sag_comp2 - sag_comp4).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Third sag to cor
    current_tr = tsearchn(tri.deform_tri_sag.X,tri.deform_tri_ax.Triangulation,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);
    b          = cartToBary(tri.control_tri, current_tr,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(tri.control_tri, current_tr, b); 
    
    b2           = cartToBary(tri.deform_tri_cor, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(tri.deform_tri_cor, current_tr, b2);     
    
    
    sag_comp3   = baryToCart(reg.target_tri_ax, current_tr, b); 
    b2          = cartToBary(tri.deform_tri_cor, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(tri.deform_tri_cor, current_tr, b2);
    
    
    error5 = sum((sag_comp2 - sag_comp4).^2);
    
    %% Fourth cor to sag
    current_tr = tsearchn(tri.deform_tri_cor.X,tri.deform_tri_cor.Triangulation,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);
    b          = cartToBary(tri.control_tri, current_tr,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(tri.control_tri, current_tr, b); 
    
    b2           = cartToBary(tri.deform_tri_sag, current_tr,[sag_comp(:,1),sag_comp(:,2),sag_comp(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(tri.deform_tri_sag, current_tr, b2);     
    
    
    sag_comp3   = baryToCart(reg.target_tri_cor, current_tr, b); 
    b2           = cartToBary(tri.deform_tri_sag, current_tr,[sag_comp3(:,1),sag_comp3(:,2),sag_comp3(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp4   = baryToCart(tri.deform_tri_sag, current_tr, b2);
    
    error6 = sum((sag_comp2 - sag_comp4).^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     current_tr = tsearchn(reg.target_tri_sag.X,reg.target_tri_sag.Triangulation,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);
%     b          = cartToBary(tri.deform_tri_ax, current_tr,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);% barycentric coordinates of the points wrt source tr
%     ax_comp   = baryToCart(tri.deform_tri_ax, current_tr, b); 
%     
%     current_tr2 = tsearchn(tri.control_tri.X,tri.control_tri.Triangulation,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);
%     b2          = cartToBary(tri.deform_tri_ax, current_tr2,[tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3)]);% barycentric coordinates of the points wrt source tr
%     ax_comp2   = baryToCart(tri.deform_tri_ax, current_tr2, b2); 
%     
%     error2 = sum((ax_comp2 - ax_comp).^2);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Third ax to cor
%     current_tr = tsearchn(reg.target_tri_ax.X,reg.target_tri_ax.Triangulation,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);
%     b          = cartToBary(tri.deform_tri_cor, current_tr,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);% barycentric coordinates of the points wrt source tr
%     cor_comp   = baryToCart(tri.deform_tri_cor, current_tr, b); 
%     
%     current_tr2 = tsearchn(tri.control_tri.X,tri.control_tri.Triangulation,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);
%     b2          = cartToBary(tri.deform_tri_cor, current_tr2,[tri.xc2_ax(:,1),tri.xc2_ax(:,2),tri.xc2_ax(:,3)]);% barycentric coordinates of the points wrt source tr
%     cor_comp2   = baryToCart(tri.deform_tri_cor, current_tr2, b2); 
%     
%     error3 = sum((cor_comp2 - cor_comp).^2);
%     
%     %% Fourth cor to ax
%     current_tr = tsearchn(reg.target_tri_cor.X,reg.target_tri_cor.Triangulation,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);
%     b          = cartToBary(tri.deform_tri_ax, current_tr,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);% barycentric coordinates of the points wrt source tr
%     ax_comp   = baryToCart(tri.deform_tri_ax, current_tr, b); 
%     
%     current_tr2 = tsearchn(tri.control_tri.X,tri.control_tri.Triangulation,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);
%     b2          = cartToBary(tri.deform_tri_ax, current_tr2,[tri.xc2_cor(:,1),tri.xc2_cor(:,2),tri.xc2_cor(:,3)]);% barycentric coordinates of the points wrt source tr
%     ax_comp2   = baryToCart(tri.deform_tri_ax, current_tr2, b2); 
%     
%     error4 = sum((ax_comp2 - ax_comp).^2);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Fifth sag to cor
%     current_tr = tsearchn(reg.target_tri_sag.X,reg.target_tri_sag.Triangulation,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);
%     b          = cartToBary(tri.deform_tri_cor, current_tr,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);% barycentric coordinates of the points wrt source tr
%     cor_comp   = baryToCart(tri.deform_tri_cor, current_tr, b); 
%     
%     current_tr2 = tsearchn(tri.control_tri.X,tri.control_tri.Triangulation,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);
%     b2          = cartToBary(tri.deform_tri_cor, current_tr2,[tri.xc3_sag(:,1),tri.xc3_sag(:,2),tri.xc3_sag(:,3)]);% barycentric coordinates of the points wrt source tr
%     cor_comp2   = baryToCart(tri.deform_tri_cor, current_tr2, b2); 
%     
%     error5 = sum((cor_comp2 - cor_comp).^2);
%     
%     %% Sixth cor to sag
%     current_tr = tsearchn(reg.target_tri_cor.X,reg.target_tri_cor.Triangulation,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);
%     b          = cartToBary(tri.deform_tri_sag, current_tr,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);% barycentric coordinates of the points wrt source tr
%     ax_comp   = baryToCart(tri.deform_tri_sag, current_tr, b); 
%     
%     current_tr2 = tsearchn(tri.control_tri.X,tri.control_tri.Triangulation,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);
%     b2          = cartToBary(tri.deform_tri_sag, current_tr2,[tri.xc3_cor(:,1),tri.xc3_cor(:,2),tri.xc3_cor(:,3)]);% barycentric coordinates of the points wrt source tr
%     ax_comp2   = baryToCart(tri.deform_tri_sag, current_tr2, b2); 
%     
%     error6 = sum((ax_comp2 - ax_comp).^2);
%     
% %     plot3(tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3),'g+');
% %     plot3(sag_comp2(:,1),sag_comp2(:,2),sag_comp2(:,3),'g*');
% %     plot3(sag_comp(:,1),sag_comp(:,2),sag_comp(:,3),'r*');
%     
%     error = ( error1 + error2 + error3 + error4 )./(6*34);
    error = sum(error1 + error2 + error3 + error4 + error5 + error6)/(2 * size(tri.xc1_ax,1) + 2* size(tri.xc2_ax,1) + 2* size(tri.xc3_sag,1));
    save(save_name,'error');
    disp(['Error: ', num2str(error)]);
end

error