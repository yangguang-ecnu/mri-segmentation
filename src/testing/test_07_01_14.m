clc;
close all;

X_ax = [];
Y_ax = [];
Z_ax = [];

X_sag = [];
Y_sag = [];
Z_sag = [];
%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1);hold on

for ax=1:size(views.axial,3)-1:size(views.axial,3)
    M = zeros(4,4);
    
    M(1:3,4) = views.axial_info{ax}.ImagePositionPatient;
    M(4,4) = 1;
    M(1:3,1) = views.axial_info{ax}.ImageOrientationPatient(1:3).*views.axial_info{ax}.PixelSpacing(1);
    M(1:3,2) = views.axial_info{ax}.ImageOrientationPatient(4:6).*views.axial_info{ax}.PixelSpacing(2);
    
    i = [0 511];
    j = [0 511];
    
    x = [];
    y = [];
    z = [];
    
    
    for k=1:length(i)
        for l=1:length(j)
            p = M*[i(k) j(l) 0 1]';
            
            x = [x p(1)];
            y = [y p(2)];
            z = [z p(3)];
            
        end
    end
    
    P1 = [x(1) y(1) z(1)];
    P2 = [x(2) y(2) z(2)];
    P3 = [x(3) y(3) z(3)];
    P4 = [x(4) y(4) z(4)];
    
    x = [P1(1) P4(1) P2(1)];   %# [xorigin xA xB] coordinates in 3-D space
    y = [P1(2) P4(3) P2(2)];   %# [yorigin yA yB] coordinates in 3-D space
    z = [P1(3) P4(3) P2(3)];   %# [zorigin zA zB] coordinates in 3-D space
    
    origin = [512 1];  %# Vertex of triangle in image space
    U = [0 511];       %# Vector from origin to point A in image space
    V = [-511 511];
    img = views.axial(:,:,ax);  %# Sample image for texture map
    
    A = origin+U;  %# Point A
    B = origin+V;  %# Point B
    C = B-U;     %# Point C
    
    [nRows,nCols,nPages] = size(img);  %# Image dimensions
    inputCorners = [origin; ...        %# Corner coordinates of input space
        A; ...
        B; ...
        C];
    
    outputCorners = [origin; ...      %# Corner coordinates of output space
        A; ...
        B; ...
        C];
    
    tform = maketform('projective',...  %# Make the transformation structure
        inputCorners,...
        outputCorners);
    
    % triTexture = imtransform(img,tform,'bicubic',...  %# Transform the image
    %                          'xdata',[1 nCols],...
    %                          'ydata',[1 nRows],...
    %                          'size',[nRows nCols]);
    
    x = [P3(1) P4(1) P2(1) P1(1)];   %# [xorigin xA xB] coordinates in 3-D space
    y = [P3(2) P4(2) P2(2) P1(2)];   %# [yorigin yA yB] coordinates in 3-D space
    z = [P3(3) P4(3) P2(3) P1(3)];   %# [zorigin zA zB] coordinates in 3-D space
    
    
    index = [4 1;3 2];  %# Index used to create 2-by-2 surface coordinates
    X1 = x(index);      %# x coordinates of surface
    Y1 = y(index);       %# y coordinates of surface
    Z1 = z(index);        %# z coordinates of surface
    
    
    points = [x' y' z'];
    fill3(points(:,1),points(:,2),points(:,3),'r');hold on
    alpha(0.5)
    
    X_ax = [X_ax x'];
    Y_ax = [Y_ax y'];
    Z_ax = [Z_ax z'];
end
%plot3([x x(1)],[y y(1)],[z z(1)]);

%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sag=1:size(views.sagittal,3)-1:size(views.sagittal,3)
    M = zeros(4,4);
    
    M(1:3,4) = views.sagittal_info{sag}.ImagePositionPatient;
    M(4,4) = 1;
    M(1:3,1) = views.sagittal_info{sag}.ImageOrientationPatient(1:3).*views.sagittal_info{sag}.PixelSpacing(1);
    M(1:3,2) = views.sagittal_info{sag}.ImageOrientationPatient(4:6).*views.sagittal_info{sag}.PixelSpacing(2);
    
    
    
    i = [0 511];
    j = [0 511];
    
    x = [];
    y = [];
    z = [];
    for k=1:2
        for l=1:2
            p = M*[j(k) i(l) 0 1]';
            
            x = [x p(1)];
            y = [y p(2)];
            z = [z p(3)];
        end
    end
    
    
    
    P1 = [x(1) y(1) z(1)];
    P2 = [x(2) y(2) z(2)];
    P3 = [x(3) y(3) z(3)];
    P4 = [x(4) y(4) z(4)];
    
    x = [P1(1) P4(1) P2(1)];   %# [xorigin xA xB] coordinates in 3-D space
    y = [P1(2) P4(3) P2(2)];   %# [yorigin yA yB] coordinates in 3-D space
    z = [P1(3) P4(3) P2(3)];   %# [zorigin zA zB] coordinates in 3-D space
    
    origin = [512 1];  %# Vertex of triangle in image space
    U = [0 511];       %# Vector from origin to point A in image space
    V = [-511 511];
    img = views.sagittal(:,:,sag);  %# Sample image for texture map
    
    A = origin+U;  %# Point A
    B = origin+V;  %# Point B
    C = B-U;     %# Point C
    
    [nRows,nCols,nPages] = size(img);  %# Image dimensions
    inputCorners = [origin; ...        %# Corner coordinates of input space
        A; ...
        B; ...
        C];
    
    outputCorners = [origin; ...      %# Corner coordinates of output space
        A; ...
        B; ...
        C];
    
    tform = maketform('projective',...  %# Make the transformation structure
        inputCorners,...
        outputCorners);
    
    x = [P3(1) P4(1) P2(1) P1(1)];   %# [xorigin xA xB] coordinates in 3-D space
    y = [P3(2) P4(2) P2(2) P1(2)];   %# [yorigin yA yB] coordinates in 3-D space
    z = [P3(3) P4(3) P2(3) P1(3)];   %# [zorigin zA zB] coordinates in 3-D space
    index = [4 1;3 2];  %# Index used to create 2-by-2 surface coordinates
    X1 = x(index);      %# x coordinates of surface
    Y1 = y(index);       %# y coordinates of surface
    Z1 = z(index);        %# z coordinates of surface
    
    points = [x' y' z'];
    fill3(points(:,1),points(:,2),points(:,3),'b');hold on
    alpha(0.3)
    
    X_sag = [X_sag x'];
    Y_sag = [Y_sag y'];
    Z_sag = [Z_sag z'];
end


%%%%%%%%%%% Plane intersection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Axial 1 Sagittal 1
N1 = cross([X_ax(1,1) Y_ax(1,1) Z_ax(1,1)]-[X_ax(2,1) Y_ax(2,1) Z_ax(2,1)],[X_ax(1,1) Y_ax(1,1) Z_ax(1,1)]-[X_ax(3,1) Y_ax(3,1) Z_ax(3,1)]);
A1 = [X_ax(4,1) Y_ax(4,1) Z_ax(4,1)];

N2 = cross([X_sag(1,1) Y_sag(1,1) Z_sag(1,1)]-[X_sag(2,1) Y_sag(2,1) Z_sag(2,1)],[X_sag(1,1) Y_sag(1,1) Z_sag(1,1)]-[X_sag(3,1) Y_sag(3,1) Z_sag(3,1)]);
A2 = [X_sag(4,1) Y_sag(4,1) Z_sag(4,1)];

[P1,dir1,check]=plane_intersect(N1,A1,N2,A2);

l = [P1-40.*dir1./norm(dir1); P1 + 40.*dir1./norm(dir1)];

plot3(l(:,1),l(:,2),l(:,3),'k');hold on
plot3(P1(1),P1(2),P1(3),'k+')

% Axial 1 Sagittal end
N1 = cross([X_ax(1,1) Y_ax(1,1) Z_ax(1,1)]-[X_ax(2,1) Y_ax(2,1) Z_ax(2,1)],[X_ax(1,1) Y_ax(1,1) Z_ax(1,1)]-[X_ax(3,1) Y_ax(3,1) Z_ax(3,1)]);
A1 = [X_ax(4,1) Y_ax(4,1) Z_ax(4,1)];

N2 = cross([X_sag(1,end) Y_sag(1,end) Z_sag(1,end)]-[X_sag(2,end) Y_sag(2,end) Z_sag(2,end)],[X_sag(1,end) Y_sag(1,end) Z_sag(1,end)]-[X_sag(3,end) Y_sag(3,end) Z_sag(3,end)]);
A2 = [X_sag(4,end) Y_sag(4,end) Z_sag(4,end)];

[P2,dir2,check]=plane_intersect(N1,A1,N2,A2);

l = [P2-40.*dir2./norm(dir2); P + 40.*dir2./norm(dir2)];

plot3(l(:,1),l(:,2),l(:,3),'k');hold on
plot3(P2(1),P2(2),P2(3),'k+')


% % Axial end Sagittal 1
N1 = cross([X_ax(1,end) Y_ax(1,end) Z_ax(1,end)]-[X_ax(2,end) Y_ax(2,end) Z_ax(2,end)],[X_ax(1,end) Y_ax(1,end) Z_ax(1,end)]-[X_ax(3,end) Y_ax(3,end) Z_ax(3,end)]);
A1 = [X_ax(4,end) Y_ax(4,end) Z_ax(4,end)];

N2 = cross([X_sag(1,1) Y_sag(1,1) Z_sag(1,1)]-[X_sag(2,1) Y_sag(2,1) Z_sag(2,1)],[X_sag(1,1) Y_sag(1,1) Z_sag(1,1)]-[X_sag(3,1) Y_sag(3,1) Z_sag(3,1)]);
A2 = [X_sag(4,1) Y_sag(4,1) Z_sag(4,1)];

[P3,dir3,check]=plane_intersect(N1,A1,N2,A2);

l = [P3-40.*dir3./norm(dir3); P3 + 40.*dir3./norm(dir3)];

plot3(l(:,1),l(:,2),l(:,3),'k');hold on
plot3(P1(1),P1(2),P1(3),'k+')

% % Axial end Sagittal end
N1 = cross([X_ax(1,end) Y_ax(1,end) Z_ax(1,end)]-[X_ax(2,end) Y_ax(2,end) Z_ax(2,end)],[X_ax(1,end) Y_ax(1,end) Z_ax(1,end)]-[X_ax(3,end) Y_ax(3,end) Z_ax(3,end)]);
A1 = [X_ax(4,end) Y_ax(4,end) Z_ax(4,end)];

N2 = cross([X_sag(1,end) Y_sag(1,end) Z_sag(1,end)]-[X_sag(2,end) Y_sag(2,end) Z_sag(2,end)],[X_sag(1,end) Y_sag(1,end) Z_sag(1,end)]-[X_sag(3,end) Y_sag(3,end) Z_sag(3,end)]);
A2 = [X_sag(4,end) Y_sag(4,end) Z_sag(4,end)];

[P4,dir4,check]=plane_intersect(N1,A1,N2,A2);

l = [P4-40.*dir4./norm(dir4); P4 + 40.*dir4./norm(dir4)];

plot3(l(:,1),l(:,2),l(:,3),'k');hold on
plot3(P1(1),P1(2),P1(3),'k+')

x = [X_ax(:,1);X_sag(:,1)];
y = [Y_ax(:,1);Y_sag(:,1)];
z = [Z_ax(:,1);Z_sag(:,1)];

V = [x y z];

[A1,b1,Aeq1,beq1] = vert2lcon([X_ax(:,1) Y_ax(:,1) Z_ax(:,1)]);
[A2,b2,Aeq2,beq2] = vert2lcon([X_sag(:,1) Y_sag(:,1) Z_sag(:,1)]);

normal = cross([X_sag(1,1) Y_sag(1,1) Z_sag(1,1)]-[X_sag(2,1) Y_sag(2,1) Z_sag(2,1)],[X_sag(1,1) Y_sag(1,1) Z_sag(1,1)]-[X_sag(3,1) Y_sag(3,1) Z_sag(3,1)]);

f = [normal -(P1(1)*normal(1) + P1(2)*normal(2) + P1(3)*normal(3))];

A = [A1;A2];
b = [b1;b2];
Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];


[x,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq);

[V,nr,nre] = qlcon2vert(x, A,b,Aeq,beq)