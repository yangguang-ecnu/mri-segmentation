function voxel2RCS(views)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Displays the 3 planes in the reference coordinate system 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Axial')
% 
% views.axial_info{1}.SliceThickness
% views.axial_info{1}.SpacingBetweenSlices
% 
% views.axial_info{1}.ImagePositionPatient
% views.axial_info{1}.ImageOrientationPatient
% 
% disp('Sagittal')
% 
% views.sagittal_info{1}.SliceThickness
% views.sagittal_info{1}.SpacingBetweenSlices
% 
% views.sagittal_info{1}.ImagePositionPatient
% views.sagittal_info{1}.ImageOrientationPatient
% 
% disp('Coronal')
% 
% views.coronal_info{1}.SliceThickness
% views.coronal_info{1}.SpacingBetweenSlices
% 
% views.coronal_info{1}.ImagePositionPatient
% views.coronal_info{1}.ImageOrientationPatient
% views.coronal_info{1}.SliceLocation


figure;
%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = zeros(4,4);

M(1:3,4) = views.axial_info{25}.ImagePositionPatient;
M(4,4) = 1;
M(1:3,1) = views.axial_info{25}.ImageOrientationPatient(1:3).*views.axial_info{25}.PixelSpacing(1);
M(1:3,2) = views.axial_info{25}.ImageOrientationPatient(4:6).*views.axial_info{25}.PixelSpacing(2);



i = [0 511];
j = [0 511];

x = [];
y = [];
z = [];
for k=1:2
    for l=1:2
    p = M*[i(k) j(l) 0 1]';
    plot3(p(1),p(2),p(3),'+r');hold on
        x = [x p(1)];
        y = [y p(2)];
        z = [z p(3)];
    end
end
plot3(x(1:2),y(1:2),z(1:2),'r');
plot3(x(1:2:3),y(1:2:3),z(1:2:3),'r');
plot3(x(2:2:4),y(2:2:4),z(2:2:4),'r');
plot3(x(3:4),y(3:4),z(3:4),'r');



%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = zeros(4,4);

M(1:3,4) = views.sagittal_info{24}.ImagePositionPatient;
M(4,4) = 1;
M(1:3,1) = views.sagittal_info{24}.ImageOrientationPatient(1:3).*views.sagittal_info{24}.PixelSpacing(1);
M(1:3,2) = views.sagittal_info{24}.ImageOrientationPatient(4:6).*views.sagittal_info{24}.PixelSpacing(2);



i = [0 511];
j = [0 511];
x = [];
y = [];
z = [];
for k=1:2
    for l=1:2
        p = M*[i(k) j(l) 0 1]';
        plot3(p(1),p(2),p(3),'+b');hold on
        x = [x p(1)];
        y = [y p(2)];
        z = [z p(3)];
    end
end
plot3(x(1:2),y(1:2),z(1:2));
plot3(x(1:2:3),y(1:2:3),z(1:2:3));
plot3(x(2:2:4),y(2:2:4),z(2:2:4));
plot3(x(3:4),y(3:4),z(3:4));
%% Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = zeros(4,4);

M(1:3,4) = views.coronal_info{25}.ImagePositionPatient;
M(4,4) = 1;
M(1:3,1) = views.coronal_info{25}.ImageOrientationPatient(1:3).*views.coronal_info{25}.PixelSpacing(1);
M(1:3,2) = views.coronal_info{25}.ImageOrientationPatient(4:6).*views.coronal_info{25}.PixelSpacing(2);



i = [0 511];
j = [0 511];
x = [];
y = [];
z = [];
for k=1:2
    for l=1:2
    p = M*[i(k) j(l) 0 1]';
    plot3(p(1),p(2),p(3),'+g');hold on
        x = [x p(1)];
        y = [y p(2)];
        z = [z p(3)];
    end
end
plot3(x(1:2),y(1:2),z(1:2),'g');
plot3(x(1:2:3),y(1:2:3),z(1:2:3),'g');
plot3(x(2:2:4),y(2:2:4),z(2:2:4),'g');
plot3(x(3:4),y(3:4),z(3:4),'g');

P1 = [x(1) y(1) z(1)];
P2 = [x(2) y(2) z(2)];
P3 = [x(3) y(3) z(3)];
P4 = [x(4) y(4) z(4)];
normal = cross(P1-P2, P1-P3);
D = -(normal(1)*P4(1) + normal(2)*P4(2)  + normal(3)*P4(3) );

[rows cols] = size(views.coronal(:,:,25));

[Y,Z] = meshgrid(1:cols, 1:rows);

X = -(1/normal(1)).*(normal(3).*Z + normal(2).*Y + D);

% figure;
% surface('XData',X, 'YData',Y-0.5, 'ZData',Z-0.5, ...
%         'CData',convert2u8(views.coronal(end:-1:1,end:-1:1,25)), 'CDataMapping','direct', ...
%         'EdgeColor','none', 'FaceColor','texturemap')
% colormap(gray)
%surf(X, Y, Z, convert2u8(views.coronal(:,:,25)), 'FaceColor', 'texturemap')