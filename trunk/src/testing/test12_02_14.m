
figure;
k_sag = 235;
ss(:,:,1) = convert2u8(lava_flex_sag(:,:,k_sag));
ss(:,:,2) = convert2u8(lava_flex_sag(:,:,k_sag));
ss(:,:,3) = convert2u8(lava_flex_sag(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(lava_sagM{k_sag}, ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

k_ax = 74;
aa(:,:,1) = convert2u8(lava_flex(:,:,k_ax));
aa(:,:,2) = convert2u8(lava_flex(:,:,k_ax));
aa(:,:,3) = convert2u8(lava_flex(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(lava_axM{k_ax}, aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

k_cor = 215;
cc(:,:,1) = convert2u8(lava_flex_cor(:,:,k_cor));
cc(:,:,2) = convert2u8(lava_flex_cor(:,:,k_cor));
cc(:,:,3) = convert2u8(lava_flex_cor(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(lava_corM{k_cor}, cc);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

%% Vol
global vol_sag_eval
global vol_ax_eval
global vol_cor_eval 
figure;
k_sag = 11;
ss(:,:,1) = convert2u8(vol_sag_eval(:,:,k_sag));
ss(:,:,2) = convert2u8(vol_sag_eval(:,:,k_sag));
ss(:,:,3) = convert2u8(vol_sag_eval(:,:,k_sag));

[X, Y, Z, triTexture] = compute_RCS_v(sag_M{k_sag}, ss);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

k_ax = 8;
aa(:,:,1) = convert2u8(vol_ax_eval(:,:,k_ax));
aa(:,:,2) = convert2u8(vol_ax_eval(:,:,k_ax));
aa(:,:,3) = convert2u8(vol_ax_eval(:,:,k_ax));

[X, Y, Z, triTexture] = compute_RCS_v(axial_M{k_ax}, aa);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

k_cor = 9;
cc(:,:,1) = convert2u8(vol_cor_eval(:,:,k_cor));
cc(:,:,2) = convert2u8(vol_cor_eval(:,:,k_cor));
cc(:,:,3) = convert2u8(vol_cor_eval(:,:,k_cor));

[X, Y, Z, triTexture] = compute_RCS_v(cor_M{k_cor}, cc);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal