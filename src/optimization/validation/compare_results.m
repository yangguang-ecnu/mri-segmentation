%% Main evaluation, compare the different points %%

clc
serie = 7;

lava_flex_n      = cellimages2mat(dcmdir.dcmPatient.Study.Series(serie,1).Images);
lava_flex_info = dcmdir.dcmPatient.Study.Series(serie,1).ImagesInfo;

[size1, size2, size3] = size(lava_flex_n);

disp('--------- Image denoising -----')
% Image denoising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size3
    lava_flex(:,:,i) = anisodiff2D(lava_flex_n(:,:,i), 20, 1/7, 30, 1);
end

%% deformation magnitude %%
sigma = [0 3 4 5 6 7 8];

t = 0:1/14:1;
root = 'deformation_';

n_x = 4;
n_y = 4;
n_z = 3;

figure;
for i=3%1:length(sigma)
    
    %% registered triangulation
    file_name = strcat(root, num2str(sigma(i)), '.mat');
    reg = load(file_name);
    
    %% original and perturbed triangulation
    st  = coordinates_M_coordinates(lava_flex_n, lava_flex, lava_flex_info);
    tri = coordinates_deformation(t, st, n_x, n_y, n_z, sigma(i));
    
    %% 6 comparisons
    
%     current_tr = tsearchn(reg.target_tri_ax.X,reg.target_tri_ax.Triangulation,[tri.var_array1(:,1) tri.var_array1(:,2) tri.var_array1(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
    current_tr = tsearchn(reg.target_tri_ax.X,reg.target_tri_ax.Triangulation,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);
    b          = cartToBary(tri.deform_tri_sag, current_tr,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp   = baryToCart(tri.deform_tri_sag, current_tr, b); 
    
    current_tr2 = tsearchn(tri.control_tri.X,tri.control_tri.Triangulation,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);
    b2          = cartToBary(tri.deform_tri_sag, current_tr2,[tri.xc1_ax(:,1),tri.xc1_ax(:,2),tri.xc1_ax(:,3)]);% barycentric coordinates of the points wrt source tr
    sag_comp2   = baryToCart(tri.deform_tri_sag, current_tr2, b2); 
    
%     plot3(tri.xc1_sag(:,1),tri.xc1_sag(:,2),tri.xc1_sag(:,3),'g+');
    plot3(sag_comp2(:,1),sag_comp2(:,2),sag_comp2(:,3),'g*');
    plot3(sag_comp(:,1),sag_comp(:,2),sag_comp(:,3),'r*');
    
    error = sum((sag_comp2 - sag_comp).^2);
end

error
