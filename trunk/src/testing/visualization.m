%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;

target_tri_ax = TriRep(source_tri.Triangulation,xfinal(1:in.size_source,:));
target_tri_sag = TriRep(source_tri.Triangulation,xfinal(in.size_source+1:end,:));

vol_ax  = views.axial;
vol_sag = views.sagittal;

%% Axial %%
k = 9;
for i = 1:size(views.axial,1)
    for j = 1:size(views.axial,2)
        
        p_3d = axial_m{k} * [j i 1]';
        
        current_tr = tsearchn(source_tri.X,source_tri.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_ax, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_ax = axial_m1{k} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1ax = floor(tmp_v1_ax(1) + 1);
            fl2ax = floor(tmp_v1_ax(2) + 1);
            cl1ax = ceil(tmp_v1_ax(1)  + 1);
            cl2ax = ceil(tmp_v1_ax(2)  + 1);
            
            min_max_r1ax = min(max(fl2ax,1),in.rows);
            min_max_r2ax = min(max(cl2ax,1),in.rows);
            min_max_c1ax = min(max(fl1ax,1),in.cols);
            min_max_c2ax = min(max(cl1ax,1),in.cols);
            
            neigax = [vol_ax(min_max_r1ax, min_max_c1ax, k) vol_ax(min_max_r1ax, min_max_c2ax, k);...
                      vol_ax(min_max_r2ax, min_max_c1ax, k) vol_ax(min_max_r2ax, min_max_c2ax, k)];
            
            new_axial(i,j,k) = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax));
        else
            new_axial(i,j,k) = 0;
        end
        
    end
end

figure;
subplot(121);imshow(vol_ax(:,:,k),   []);title('Original');
subplot(122);imshow(new_axial(:,:,k),[]);title('After Deformation');
%% Sagittal %%

k = 9;
for i = 1:size(views.sagittal,1)
    for j = 1:size(views.sagittal,2)
        
        p_3d = sag_m{k} * [j i 1]';
        
        current_tr = tsearchn(source_tri.X,source_tri.Triangulation,[p_3d(1) p_3d(2) p_3d(3)]); % calculate the tetrahedron where p_3d belongs
        
        if ~isnan(current_tr)
            c2b_coord = cartToBary(source_tri,current_tr,[p_3d(1) p_3d(2) p_3d(3)]); % get the barycentric coordinates
            b2c_ncoord = baryToCart(target_tri_sag, current_tr, c2b_coord); % get the cartesian coordinates
            
            tmp_v1_sag = sag_m1{k} * [b2c_ncoord 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            fl1sag = floor(tmp_v1_sag(1) + 1);
            fl2sag = floor(tmp_v1_sag(2) + 1);
            cl1sag = ceil(tmp_v1_sag(1)  + 1);
            cl2sag = ceil(tmp_v1_sag(2)  + 1);
            
            min_max_r1sag = min(max(fl2sag,1),in.rows);
            min_max_r2sag = min(max(cl2sag,1),in.rows);
            min_max_c1sag = min(max(fl1sag,1),in.cols);
            min_max_c2sag = min(max(cl1sag,1),in.cols);
            
            neigsag = [vol_sag(min_max_r1sag, min_max_c1sag, k) vol_sag(min_max_r1sag, min_max_c2sag, k);...
                       vol_sag(min_max_r2sag, min_max_c1sag, k) vol_sag(min_max_r2sag, min_max_c2sag, k)];
            
            new_sagittal(i,j,k) = bilinear_interpolation(tmp_v1_sag(2), tmp_v1_sag(1), double(neigsag));
        else
            new_sagittal(i,j,k) = 0;
        end
        
    end
end

figure;
subplot(121);imshow(vol_sag(:,:,k),   []);title('Original');
subplot(122);imshow(new_sagittal(:,:,k),[]);title('After Deformation');


%% Visualize in the RCS
% Original
figure;
aa(:,:,1) = convert2u8(vol_ax(:,:,k));
aa(:,:,2) = convert2u8(vol_ax(:,:,k));
aa(:,:,3) = convert2u8(vol_ax(:,:,k));

[X, Y, Z, triTexture] = compute_RCS_v(axial_m{k},aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag(:,:,k));
ss(:,:,1) = convert2u8(vol_sag(:,:,k));
ss(:,:,3) = convert2u8(vol_sag(:,:,k));

[X, Y, Z, triTexture] = compute_RCS_v(sag_m{k},ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal
% Deformation
figure;
aa(:,:,1) = convert2u8(new_axial(:,:,k));
aa(:,:,2) = convert2u8(new_axial(:,:,k));
aa(:,:,3) = convert2u8(new_axial(:,:,k));

[X, Y, Z, triTexture] = compute_RCS_v(axial_m{k},aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

ss(:,:,2) = convert2u8(new_sagittal(:,:,k));
ss(:,:,1) = convert2u8(new_sagittal(:,:,k));
ss(:,:,3) = convert2u8(new_sagittal(:,:,k));

[X, Y, Z, triTexture] = compute_RCS_v(sag_m{k},ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal
% TO DO:
%% Coronal %%