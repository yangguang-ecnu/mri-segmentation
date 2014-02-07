total_ax = size(views.axial,3);
total_sag = size(views.sagittal,3);

axial_m  = cell(1,total_ax);
axial_m1 = cell(1,total_ax);

X_ax = [];Y_ax = [];Z_ax = [];
X_sag = [];Y_sag = [];Z_sag = [];

vol_ax  = views.axial;
vol_sag = views.sagittal;

for ax=1:size(views.axial,3)
    
    [axial_m{ax}, axial_m1{ax}, tr_ax] = compute_M_M1_ortho(views.axial_info{ax}, 0 ,1);
    
    [x,y,z] = calculate4cornersk( axial_m{ax},ax );

    X_ax = [X_ax x'];
    Y_ax = [Y_ax y'];
    Z_ax = [Z_ax z'];

    
    if ax == 1
        N1 = cross([X_ax(1,ax) Y_ax(1,ax) Z_ax(1,ax)]-[X_ax(2,ax) Y_ax(2,ax) Z_ax(2,ax)],[X_ax(1,ax) Y_ax(1,ax) Z_ax(1,ax)]-[X_ax(3,ax) Y_ax(3,ax) Z_ax(3,ax)]); % normal to the axial (ax)

        N1 = N1./norm(N1)

        N1 = [0 0 1];
        axial_m{ax}
    end

%     plane_ax(ax,:) = [N1 -(N1(1)*X_ax(4,ax) + N1(2)*Y_ax(4,ax) + N1(3)*Z_ax(4,ax))]; % (A,B,C,D) of the plane Ax + By + Cz + D =0
    
end
disp('--------- Calculate M^(-1) and plane eq. for each slice in the second direction -----')
%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sag_m = cell(1,total_sag);
sag_m1 = cell(1,total_sag);

for sag=1:size(views.sagittal,3)
    
    [sag_m{sag}, sag_m1{sag}, tr_sag] = compute_M_M1_ortho(views.sagittal_info{sag}, 0, 2);

    [x,y,z] = calculate4cornersk( sag_m{sag} ,sag);
       
    X_sag = [X_sag x'];
    Y_sag = [Y_sag y'];
    Z_sag = [Z_sag z'];
    
    if sag == 1
        N2 = cross([X_sag(1,sag) Y_sag(1,sag) Z_sag(1,sag)]-[X_sag(2,sag) Y_sag(2,sag) Z_sag(2,sag)],[X_sag(1,sag) Y_sag(1,sag) Z_sag(1,sag)]-[X_sag(3,sag) Y_sag(3,sag) Z_sag(3,sag)]);
        N2 = N2./norm(N2)

        N2 = [-1 0 0];
        sag_m{sag}
    end

%     plane_sag(sag,:) = [N2 -(N2(1)*X_sag(4,sag) + N2(2)*Y_sag(4,sag) + N2(3)*Z_sag(4,sag))]; % (A,B,C,D) of the plane Ax + By + Cz + D =0
    
end

figure;
fill3(X_ax(:,1),Y_ax(:,1),Z_ax(:,1),'r');hold on % first axial plane
fill3(X_ax(:,total_ax),Y_ax(:,total_ax),Z_ax(:,total_ax),'r');hold on % first axial plane
fill3(X_sag(:,1),Y_sag(:,1),Z_sag(:,1),'b');hold on % first sagittal plane
fill3(X_sag(:,total_sag),Y_sag(:,total_sag),Z_sag(:,total_sag),'b');hold on % second sagittal plane
alpha(.2)

k = 10;

figure;
aa(:,:,1) = convert2u8(vol_ax(:,:,k));
aa(:,:,2) = convert2u8(vol_ax(:,:,k));
aa(:,:,3) = convert2u8(vol_ax(:,:,k));

[X, Y, Z, triTexture] = compute_RCS_v2(axial_m{k},aa,k);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
               'FaceColor','texturemap',...
               'EdgeColor','none');hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag(:,:,k));
ss(:,:,1) = convert2u8(vol_sag(:,:,k));
ss(:,:,3) = convert2u8(vol_sag(:,:,k));

[X, Y, Z, triTexture] = compute_RCS_v2(sag_m{k},ss,k);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag(:,:,end));
ss(:,:,1) = convert2u8(vol_sag(:,:,end));
ss(:,:,3) = convert2u8(vol_sag(:,:,end));

[X, Y, Z, triTexture] = compute_RCS_v2(sag_m{end},ss,25);
hSurface = surf(X,Y,Z,triTexture,...          %# Plot texture-mapped surface
              'FaceColor','texturemap',...
              'EdgeColor','none');hold on
axis equal 