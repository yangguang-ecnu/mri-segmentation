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

for ax=1:size(views.axial,3)
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
%     fill3(points(:,1),points(:,2),points(:,3),'r');hold on
%     alpha(0.5)
    
    X_ax = [X_ax x'];
    Y_ax = [Y_ax y'];
    Z_ax = [Z_ax z'];
end


%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sag=1:size(views.sagittal,3)
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
%     fill3(points(:,1),points(:,2),points(:,3),'b');hold on
%     alpha(0.3)
    
    X_sag = [X_sag x'];
    Y_sag = [Y_sag y'];
    Z_sag = [Z_sag z'];
end


%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [];
ax=1:size(views.axial,3);
sag=1:size(views.sagittal,3);

options = optimset('Display','off');

for i=1:length(ax)
    for j=1:length(sag)
    
    
        x_t = [X_ax(:,i);X_sag(:,j)];
        y_t = [Y_ax(:,i);Y_sag(:,j)];
        z_t = [Z_ax(:,i);Z_sag(:,j)];

        [A1,b1,Aeq1,beq1] = vert2lcon([X_ax(:,i) Y_ax(:,i) Z_ax(:,i)]);
        [A2,b2,Aeq2,beq2] = vert2lcon([X_sag(:,j) Y_sag(:,j) Z_sag(:,j)]);

        normal = cross([X_sag(1,j) Y_sag(1,j) Z_sag(1,j)]-[X_sag(2,j) Y_sag(2,j) Z_sag(2,j)],[X_sag(1,j) Y_sag(1,j) Z_sag(1,j)]-[X_sag(3,j) Y_sag(3,j) Z_sag(3,j)]);

        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];


        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options);

        [V1,nr,nre] = qlcon2vert(x0, A,b,Aeq,beq);

        %plot3(V1(:,1),V1(:,2),V1(:,3),'g*');hold on

        %a = line(V1(:,1),V1(:,2),V1(:,3));
        %line([V1(1,1) V1(1,2) V1(1,3)], [V1(2,1),V1(2,2),V1(2,3)], 'g')

        t = 0:.5:1;
        vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];

        for k=1:length(t)
            plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'k+');hold on
            var = [var;V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)];
            vert_mat(k,j,i) = (i-1)*(length(ax)-1) + (k-1)+ j-1 + (k-1)*(length(sag)-1);
        end
        
        
    end
end

%[dt] = MyCrust(var);
dt = DelaunayTri(var);%delaunay(var);%
FV2.faces = dt;
FV2.vertices = var;
trisurf(FV2.faces,FV2.vertices(:,1),FV2.vertices(:,2),FV2.vertices(:,3),'facecolor','r','facealpha',.1,'edgecolor','b','facelighting','flat');camlight
axis equal

% edgeIndex = edges(dt);              %# Triangulation edge indices
% midpts = [mean(var(edgeIndex,1),2) ...  %# Triangulation edge midpoints
%           mean(var(edgeIndex,2),2) ...
%           mean(var(edgeIndex,3),2)];
% nearIndex = nearestNeighbor(dt,midpts);  %# Find the vertex nearest the midpoints
% keepIndex = (nearIndex == edgeIndex(:,1)) | ...  %# Find the edges where the
%             (nearIndex == edgeIndex(:,2));                                    %#   another vertex than it is
%                                                  %#   to one of its end vertices
% edgeIndex = edgeIndex(keepIndex,:);      %# The "good" edges
% And now edgeIndex is an N-by-2 matrix where each row contains the indices into x and y for one edge that defines a "near" connection. The following plot illustrates the Delaunay triangulation (red lines), Voronoi diagram (blue lines), midpoints of the triangulation edges (black asterisks), and the "good" edges that remain in edgeIndex (thick red lines):
% 
% triplot(dt,'r');  %# Plot the Delaunay triangulation
% hold on;          %# Add to the plot
% plot3(var(edgeIndex,1).',var(edgeIndex,2).',var(edgeIndex,3).','r-','LineWidth',3);  %# Plot the "good" edges

%tetramesh(dt);