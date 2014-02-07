%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;

target_tri_ax  = TriRep(source_tri.Triangulation, xfinal(1:in.size_source,:));
target_tri_sag = TriRep(source_tri.Triangulation, xfinal(in.size_source+1:2*in.size_source,:));
target_tri_cor = TriRep(source_tri.Triangulation, xfinal(2*in.size_source+1:end,:));

k_ax  = 8;
k_sag = 9;
k_cor = 9;

%% Axial %%

for i = 1:size(views.axial,1)
    for j = 1:size(views.axial,2)
        
        p_3d = axial_m{k_ax} * [j i 1]';
        
        current_tr = tsearchn(source_tri.X,source_tri.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_ax, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_ax = axial_m1{k_ax} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1ax = floor(tmp_v1_ax(1) + 1);
            fl2ax = floor(tmp_v1_ax(2) + 1);
            cl1ax = ceil(tmp_v1_ax(1)  + 1);
            cl2ax = ceil(tmp_v1_ax(2)  + 1);
            
            min_max_r1ax = min(max(fl2ax,1),in.rows);
            min_max_r2ax = min(max(cl2ax,1),in.rows);
            min_max_c1ax = min(max(fl1ax,1),in.cols);
            min_max_c2ax = min(max(cl1ax,1),in.cols);
            
            neigax = [vol_ax(min_max_r1ax, min_max_c1ax, k_ax) vol_ax(min_max_r1ax, min_max_c2ax, k_ax);...
                      vol_ax(min_max_r2ax, min_max_c1ax, k_ax) vol_ax(min_max_r2ax, min_max_c2ax, k_ax)];
            
            new_axial(i,j,k_ax) = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax));
        else
            new_axial(i,j,k_ax) = views.axial(i,j,k_ax);
        end
        
    end
end

% figure;
% subplot(121);imshow(vol_ax(:,:,k),   []);title('Original');
% subplot(122);imshow(new_axial(:,:,k),[]);title('After Deformation');
%% Sagittal %%

for i = 1:size(views.sagittal,1)
    for j = 1:size(views.sagittal,2)
        
        p_3d = sag_m{k_sag} * [j i 1]';
        
        current_tr = tsearchn(source_tri.X,source_tri.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_sag, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_sag = sag_m1{k_sag} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1sag = floor(tmp_v1_sag(1) + 1);
            fl2sag = floor(tmp_v1_sag(2) + 1);
            cl1sag = ceil(tmp_v1_sag(1)  + 1);
            cl2sag = ceil(tmp_v1_sag(2)  + 1);
            
            min_max_r1sag = min(max(fl2sag,1),in.rows);
            min_max_r2sag = min(max(cl2sag,1),in.rows);
            min_max_c1sag = min(max(fl1sag,1),in.cols);
            min_max_c2sag = min(max(cl1sag,1),in.cols);
            
            neigsag = [vol_sag(min_max_r1sag, min_max_c1sag, k_sag) vol_sag(min_max_r1sag, min_max_c2sag, k_sag);...
                       vol_sag(min_max_r2sag, min_max_c1sag, k_sag) vol_sag(min_max_r2sag, min_max_c2sag, k_sag)];
            
            new_sagittal(i,j,k_sag) = bilinear_interpolation(tmp_v1_sag(2), tmp_v1_sag(1), double(neigsag));
        else
            new_sagittal(i,j,k_sag) = views.sagittal(i,j,k_sag);
        end
        
    end
end

% figure;
% subplot(121);imshow(vol_sag(:,:,k),   []);title('Original');
% subplot(122);imshow(new_sagittal(:,:,k),[]);title('After Deformation');

%% Coronal %%

for i = 1:size(views.coronal,1)
    for j = 1:size(views.coronal,2)
        
        p_3d = cor_m{k_cor} * [j i 1]';
        
        current_tr = tsearchn(source_tri.X,source_tri.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_cor, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_cor = cor_m1{k_cor} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1cor = floor(tmp_v1_cor(1) + 1);
            fl2cor = floor(tmp_v1_cor(2) + 1);
            cl1cor = ceil(tmp_v1_cor(1)  + 1);
            cl2cor = ceil(tmp_v1_cor(2)  + 1);
            
            min_max_r1cor = min(max(fl2cor,1),in.rows);
            min_max_r2cor = min(max(cl2cor,1),in.rows);
            min_max_c1cor = min(max(fl1cor,1),in.cols);
            min_max_c2cor = min(max(cl1cor,1),in.cols);
            
            neigcor = [vol_cor(min_max_r1cor, min_max_c1cor, k_cor) vol_cor(min_max_r1cor, min_max_c2cor, k_cor);...
                       vol_cor(min_max_r2cor, min_max_c1cor, k_cor) vol_cor(min_max_r2cor, min_max_c2cor, k_cor)];
            
            new_coronal(i,j,k_cor) = bilinear_interpolation(tmp_v1_cor(2), tmp_v1_cor(1), double(neigcor));
        else
            new_coronal(i,j,k_cor) = views.coronal(i,j,k_cor);
        end
        
    end
end

% figure;
% subplot(121);imshow(vol_cor(:,:,k),   []);title('Original');
% subplot(122);imshow(new_coronal(:,:,k),[]);title('After Deformation');

%% Visualize in the RCS
% Original Axial and Sagittal
figure;
aa(:,:,1) = convert2u8(vol_ax(:,:,k_ax));
aa(:,:,2) = convert2u8(vol_ax(:,:,k_ax));
aa(:,:,3) = convert2u8(vol_ax(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_m{k_ax},aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag(:,:,k_sag));
ss(:,:,1) = convert2u8(vol_sag(:,:,k_sag));
ss(:,:,3) = convert2u8(vol_sag(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_m{k_sag},ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal

% Original Coronal and Sagittal
figure;
% k_cor = 1;
cc(:,:,2) = convert2u8(vol_cor(:,:,k_cor));
cc(:,:,1) = convert2u8(vol_cor(:,:,k_cor));
cc(:,:,3) = convert2u8(vol_cor(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_m{k_cor},cc);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal
% k_sag = 24;
ss(:,:,2) = convert2u8(vol_sag(:,:,k_sag));
ss(:,:,1) = convert2u8(vol_sag(:,:,k_sag));
ss(:,:,3) = convert2u8(vol_sag(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_m{k_sag},ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal  

% Original Coronal and Axial
figure;
aa(:,:,1) = convert2u8(vol_ax(:,:,k_ax));
aa(:,:,2) = convert2u8(vol_ax(:,:,k_ax));
aa(:,:,3) = convert2u8(vol_ax(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_m{k_ax},aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

cc(:,:,2) = convert2u8(vol_cor(:,:,k_cor));
cc(:,:,1) = convert2u8(vol_cor(:,:,k_cor));
cc(:,:,3) = convert2u8(vol_cor(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_m{k_cor},cc);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal


% Deformation Axial and Sagittal
figure;
aa(:,:,1) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,2) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,3) = convert2u8(new_axial(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_m{k_ax},aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

ss(:,:,2) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,1) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,3) = convert2u8(new_sagittal(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_m{k_sag},ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal

% Deformation Axial and Coronal
figure;

aa(:,:,1) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,2) = convert2u8(new_axial(:,:,k_ax));
aa(:,:,3) = convert2u8(new_axial(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_m{k_ax},aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

cc(:,:,2) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,1) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,3) = convert2u8(new_coronal(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_m{k_cor},cc);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal

% Deformation Sagittal and Coronal
figure;

cc(:,:,2) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,1) = convert2u8(new_coronal(:,:,k_cor));
cc(:,:,3) = convert2u8(new_coronal(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_m{k_cor},cc);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal


ss(:,:,2) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,1) = convert2u8(new_sagittal(:,:,k_sag));
ss(:,:,3) = convert2u8(new_sagittal(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_m{k_sag},ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal
