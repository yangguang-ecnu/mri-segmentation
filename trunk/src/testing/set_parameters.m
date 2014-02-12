clc

set_up = 1;

X_ax = [];
Y_ax = [];
Z_ax = [];

X_sag = [];
Y_sag = [];
Z_sag = [];

X_cor = [];
Y_cor = [];
Z_cor = [];


%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global vol_ax
global vol_sag
global vol_cor

global new_im_ax
global new_im_sag
global new_im_cor

global gradx_ax
global grady_ax
global gradx_sag
global grady_sag
global gradx_cor
global grady_cor

vol_ax  = double(views.axial);
vol_sag = double(views.sagittal);
vol_cor = double(views.coronal);

new_im_ax = zeros(size(views.axial));
new_im_sag = zeros(size(views.sagittal));

rows = size(views.axial,1);
cols = size(views.axial,2);

total_ax = size(views.axial,3);
total_sag = size(views.sagittal,3);
total_cor = size(views.coronal,3);

global plane_ax
global plane_sag
global plane_cor

plane_ax = zeros(total_ax,4);
plane_sag = zeros(total_sag,4);
plane_cor = zeros(total_cor,4);

global axial_m
global axial_m1

axial_m = cell(1,total_ax);
axial_m1 = cell(1,total_ax);

global t
t = 0:1/4:1;

global ortho
ortho = 1;

options = optimset('Display','off');

disp('--------- Image denoising -----')
%% Image denoising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:total_ax
    vol_ax(:,:,i) = anisodiff2D(views.axial(:,:,i), 10, 1/7, 30, 1);
end

for i=1:total_sag
    vol_sag(:,:,i) = anisodiff2D(views.sagittal(:,:,i), 10, 1/7, 30, 1);
end
for i=1:total_cor
    vol_cor(:,:,i) = anisodiff2D(views.coronal(:,:,i), 10, 1/7, 30, 1);
end


disp('--------- Calculate the image gradients -----')
%% Image gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:total_ax
    [gradx_ax(:,:,i)  grady_ax(:,:,i)]  = gradient(vol_ax(:,:,i));
end

for i=1:total_sag
    [gradx_sag(:,:,i) grady_sag(:,:,i)] = gradient(vol_sag(:,:,i));
end
for i=1:total_cor
    [gradx_cor(:,:,i) grady_cor(:,:,i)] = gradient(vol_cor(:,:,i));
end

disp('--------- Calculate M^(-1) and plane eq. for each slice in the first direction -----')
M = zeros(4,3);

for ax=1:size(views.axial,3)
    
    [axial_m{ax}, axial_m1{ax}, tr_ax] = compute_M_M1(views.axial_info{ax}, ortho, 1);
    
    [x,y,z] = calculate4corners( axial_m{ax} );

    X_ax = [X_ax x'];
    Y_ax = [Y_ax y'];
    Z_ax = [Z_ax z'];
     
    if ax == 1
        N1 = cross([X_ax(1,ax) Y_ax(1,ax) Z_ax(1,ax)]-[X_ax(2,ax) Y_ax(2,ax) Z_ax(2,ax)],[X_ax(1,ax) Y_ax(1,ax) Z_ax(1,ax)]-[X_ax(3,ax) Y_ax(3,ax) Z_ax(3,ax)]); % normal to the axial (ax)
        N1 = N1./norm(N1)
    end

    plane_ax(ax,:) = [N1 -(N1(1)*X_ax(4,ax) + N1(2)*Y_ax(4,ax) + N1(3)*Z_ax(4,ax))]; % (A,B,C,D) of the plane Ax + By + Cz + D =0
    
end
disp('--------- Calculate M^(-1) and plane eq. for each slice in the second direction -----')
%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global sag_m
global sag_m1

sag_m = cell(1,total_sag);
sag_m1 = cell(1,total_sag);

for sag=1:size(views.sagittal,3)
    
    [sag_m{sag}, sag_m1{sag}, tr_sag] = compute_M_M1(views.sagittal_info{sag}, ortho, 2);

    [x,y,z] = calculate4corners( sag_m{sag} );
       
    X_sag = [X_sag x'];
    Y_sag = [Y_sag y'];
    Z_sag = [Z_sag z'];
    
    if sag == 1
        N2 = cross([X_sag(1,sag) Y_sag(1,sag) Z_sag(1,sag)]-[X_sag(2,sag) Y_sag(2,sag) Z_sag(2,sag)],[X_sag(1,sag) Y_sag(1,sag) Z_sag(1,sag)]-[X_sag(3,sag) Y_sag(3,sag) Z_sag(3,sag)]);
        N2 = N2./norm(N2)
    end

    plane_sag(sag,:) = [N2 -(N2(1)*X_sag(4,sag) + N2(2)*Y_sag(4,sag) + N2(3)*Z_sag(4,sag))]; % (A,B,C,D) of the plane Ax + By + Cz + D =0
    
end

disp('--------- Calculate M^(-1) and plane eq. for each slice in the third direction -----')
%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global cor_m
global cor_m1

cor_m  = cell(1,total_cor);
cor_m1 = cell(1,total_cor);

for cor=1:size(views.coronal,3)
    
    [cor_m{cor}, cor_m1{cor}] = compute_M_M1(views.coronal_info{cor}, ortho, 3);
    
    [x,y,z] = calculate4corners( cor_m{cor} );
    
    %     points = [x' y' z'];
    %     fill3(points(:,1),points(:,2),points(:,3),'b');hold on
    
    X_cor = [X_cor x'];
    Y_cor = [Y_cor y'];
    Z_cor = [Z_cor z'];
    
    if cor == 1
        N3 = cross(-[X_cor(1,cor) Y_cor(1,cor) Z_cor(1,cor)]+[X_cor(2,cor) Y_cor(2,cor) Z_cor(2,cor)],-[X_cor(1,cor) Y_cor(1,cor) Z_cor(1,cor)]+[X_cor(3,cor) Y_cor(3,cor) Z_cor(3,cor)]);
        N3 = N3./norm(N3) 
    end
    
    plane_cor(cor,:) = [N3 -(N3(1)*X_cor(4,cor) + N3(2)*Y_cor(4,cor) + N3(3)*Z_cor(4,cor))]; % (A,B,C,D) of the plane Ax + By + Cz + D =0
    
end

figure;
fill3(X_ax(:,1),Y_ax(:,1),Z_ax(:,1),'r');hold on % first axial plane
fill3(X_ax(:,total_ax),Y_ax(:,total_ax),Z_ax(:,total_ax),'r');hold on % first axial plane
fill3(X_sag(:,1),Y_sag(:,1),Z_sag(:,1),'b');hold on % first sagittal plane
fill3(X_sag(:,total_sag),Y_sag(:,total_sag),Z_sag(:,total_sag),'b');hold on % second sagittal plane
alpha(.2)

figure;
fill3(X_ax(:,1),Y_ax(:,1),Z_ax(:,1),'r');hold on % first axial plane
fill3(X_ax(:,total_ax),Y_ax(:,total_ax),Z_ax(:,total_ax),'r');hold on % first axial plane
fill3(X_cor(:,1),Y_cor(:,1),Z_cor(:,1),'b');hold on % first sagittal plane
fill3(X_cor(:,total_cor),Y_cor(:,total_cor),Z_cor(:,total_cor),'b');hold on % second sagittal plane
alpha(.2)

figure;
fill3(X_cor(:,1),Y_cor(:,1),Z_cor(:,1),'r');hold on % first sagittal plane
fill3(X_cor(:,total_cor),Y_cor(:,total_cor),Z_cor(:,total_cor),'r');hold on % second sagittal plane
fill3(X_sag(:,1),Y_sag(:,1),Z_sag(:,1),'b');hold on % first sagittal plane
fill3(X_sag(:,total_sag),Y_sag(:,total_sag),Z_sag(:,total_sag),'b');hold on % second sagittal plane
alpha(.2)
disp('--------- Calculate the s1 x s2 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax  = 1:size(views.axial,3);
sag = 1:size(views.sagittal,3);
cor = 1:size(views.coronal,3);

global var_cell1
global var_array1

var_cell1  = cell(size(views.axial,3),size(views.sagittal,3),length(t));
var_array1 = zeros(size(views.axial,3)*size(views.sagittal,3)*length(t),3);

for i=1:length(ax)
    
    for j=1:length(sag)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([X_ax(:,i) Y_ax(:,i) Z_ax(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([X_sag(:,j) Y_sag(:,j) Z_sag(:,j)]); %% constraints of the sagittal plane
        
        normal = cross([X_sag(1,j) Y_sag(1,j) Z_sag(1,j)]-[X_sag(2,j) Y_sag(2,j) Z_sag(2,j)],[X_sag(1,j) Y_sag(1,j) Z_sag(1,j)]-[X_sag(3,j) Y_sag(3,j) Z_sag(3,j)]); % direction of the intersection line
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        
        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
        
        % plot3(V1(:,1),V1(:,2),V1(:,3),'g*');hold on
        
        vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
        
        for k=1:length(t)
            
            % plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');%,'MarkerSize',i);hold on
            var_cell1{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
            
            ind_tmp = sub2ind([size(var_cell1,3) size(var_cell1,2) size(var_cell1,1)],k,j,i);
            %ind_tmp1 = k + length(t)*(j-1 + length(sag)*(i-1));% it is faster than sub2ind
            
            var_array1(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
            var_array1(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
            var_array1(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
            
            %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [i1, j1, i2, j2, real_v]  = compute_coord(axial_m1{i}, [var_cell1{i,j,k} 1], rows, cols);
            
            neig = [vol_ax(i1, j1, i)   vol_ax(i1, j2, i);...
                    vol_ax(i2, j1, i)   vol_ax(i2, j2, i)];
            
            new_im_ax(i1, j1, i) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
            diff1 = new_im_ax(i1, j1, i);
            %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [i1, j1, i2, j2, real_v]  = compute_coord(sag_m1{j}, [var_cell1{i,j,k} 1], rows, cols);
            
            neig = [vol_sag(i1, j1, j)   vol_sag(i1, j2, j);...
                    vol_sag(i2, j1, j)   vol_sag(i2, j2, j)];
            
            new_im_sag(i1, j1, j) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
            diff2 = new_im_sag(i1, j1, j);
            
            diff(ind_tmp) = abs(diff1 - diff2);
        end
        
        
    end
end

disp('--------- Calculate the s1 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global t2
t2 = 0:1/5:1;


global var_cell2
global var_array2

var_cell2 = cell(size(views.axial,3),size(views.coronal,3),length(t));
var_array2 = zeros(size(views.axial,3)*size(views.coronal,3)*length(t),3);

options = optimset('Display','off');

figure;
for i=1:length(ax)
    
    for j=1:length(cor)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([X_ax(:,i) Y_ax(:,i) Z_ax(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([X_cor(:,j) Y_cor(:,j) Z_cor(:,j)]); %% constraints of the sagittal plane
        
        normal = cross([X_cor(1,j) Y_cor(1,j) Z_cor(1,j)]-[X_cor(2,j) Y_cor(2,j) Z_cor(2,j)],[X_cor(1,j) Y_cor(1,j) Z_cor(1,j)]-[X_cor(3,j) Y_cor(3,j) Z_cor(3,j)]); % direction of the intersection line
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        
        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
        
%         plot3(V1(:,1),V1(:,2),V1(:,3),'g*');hold on
     
        vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
        
        for k=1:length(t)
            
            plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');%,'MarkerSize',i);hold on
            var_cell2{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
            
            ind_tmp = sub2ind([size(var_cell2,3) size(var_cell2,2) size(var_cell2,1)],k,j,i);
            %ind_tmp1 = k + length(t)*(j-1 + length(sag)*(i-1));% it is faster than sub2ind
            
            var_array2(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
            var_array2(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
            var_array2(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
            
        end
        
        
    end
end


fill3(X_ax(:,1),Y_ax(:,1),Z_ax(:,1),'r');hold on % first axial plane
fill3(X_ax(:,total_ax),Y_ax(:,total_ax),Z_ax(:,total_ax),'r');hold on % first axial plane
fill3(X_cor(:,1),Y_cor(:,1),Z_cor(:,1),'b');hold on % first sagittal plane
fill3(X_cor(:,total_cor),Y_cor(:,total_cor),Z_cor(:,total_cor),'b');hold on % second sagittal plane
alpha(.2)
disp('--------- Calculate the s2 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global t3
t3 = 0:1/5:1;


global var_cell3
global var_array3

var_cell3 = cell(size(views.coronal,3),size(views.sagittal,3),length(t));
var_array3 = zeros(size(views.coronal,3)*size(views.sagittal,3)*length(t),3);

options = optimset('Display','off');
figure;
for i=1:length(cor)
    
    for j=1:length(sag)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([X_cor(:,i) Y_cor(:,i) Z_cor(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([X_sag(:,j) Y_sag(:,j) Z_sag(:,j)]); %% constraints of the sagittal plane
        
        if i==1 && j==1
            normal = cross([X_sag(1,j) Y_sag(1,j) Z_sag(1,j)]-[X_sag(2,j) Y_sag(2,j) Z_sag(2,j)],[X_sag(1,j) Y_sag(1,j) Z_sag(1,j)]-[X_sag(3,j) Y_sag(3,j) Z_sag(3,j)]); % direction of the intersection line
        end
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        
        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
        
        plot3(V1(:,1),V1(:,2),V1(:,3),'g*');hold on
        
        % a = line(V1(:,1),V1(:,2),V1(:,3));
        % line([V1(1,1) V1(1,2) V1(1,3)], [V1(2,1),V1(2,2),V1(2,3)], 'g')
        
        
        vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
        
        for k=1:length(t)
            
%             plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');
            var_cell3{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
            
            ind_tmp = sub2ind([size(var_cell3,3) size(var_cell3,2) size(var_cell3,1)],k,j,i);
            
            
            var_array3(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
            var_array3(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
            var_array3(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
            
        end
        
        
    end
end
fill3(X_cor(:,1),Y_cor(:,1),Z_cor(:,1),'r');hold on % first sagittal plane
fill3(X_cor(:,total_cor),Y_cor(:,total_cor),Z_cor(:,total_cor),'r');hold on % second sagittal plane
fill3(X_sag(:,1),Y_sag(:,1),Z_sag(:,1),'b');hold on % first sagittal plane
fill3(X_sag(:,total_sag),Y_sag(:,total_sag),Z_sag(:,total_sag),'b');hold on % second sagittal plane
alpha(.2)

tic
n_points = 1:length(t);%*size(vol_sag,3);


%% Prepare the constraints
cons_Aeq(1:length(n_points),:) = repmat(plane_ax(1,1:3),length(n_points),1);
cons_beq(1:length(n_points),:) = repmat(plane_ax(1,4),length(n_points),1);

cons_Aeq(1+length(n_points):2*length(n_points),:) = repmat(plane_sag(1,1:3),length(n_points),1);
cons_beq(1+length(n_points):2*length(n_points),:) = repmat(plane_sag(1,4), length(n_points),1);

%% Plot the points

% figure(1);
% for i = 1:length(n_points)
%     plot3(sol(i,1),sol(i,2),sol(i,3),'r*');hold on
%     plot3(var_array(n_points(i),1),var_array(n_points(i),2),var_array(n_points(i),3),'g*');hold on
%     line([sol(i,1),var_array(n_points(i),1)],[sol(i,2),var_array(n_points(i),2)],[sol(i,3),var_array(n_points(i),3)]);hold on
% end

%% Plot the planes
figure;
fill3(X_ax(:,1),Y_ax(:,1),Z_ax(:,1),'r');hold on % first axial plane
fill3(X_ax(:,total_ax),Y_ax(:,total_ax),Z_ax(:,total_ax),'r');hold on % first axial plane
fill3(X_sag(:,1),Y_sag(:,1),Z_sag(:,1),'b');hold on % first sagittal plane
fill3(X_sag(:,total_sag),Y_sag(:,total_sag),Z_sag(:,total_sag),'b');hold on % second sagittal plane
alpha(.2)
time = toc;

disp('--------- Calculate the source control points ( # N^3 ) -----')
%% Calculate the bounding box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bb = [xmin xmax;ymin ymax;zmin zmax]
global var_array

var_array = [var_array1;var_array2;var_array3];

bb = [min(min(var_array(:,1))) max(max(var_array(:,1))); ...
    min(min(var_array(:,2))) max(max(var_array(:,2))); ...
    min(min(var_array(:,3))) max(max(var_array(:,3)))];

% Create the source control points
nx = 3;
ny = 3;
nz = 3;

l_x = linspace(bb(1,1),bb(1,2),nx);
l_y = linspace(bb(2,1),bb(2,2),ny);
l_z = linspace(bb(3,1),bb(3,2),nz);

global source_control
source_control = zeros(nx * ny * nz,3);

for i = 1:nx
    for j = 1:ny
        
        tmp =  1:nz;
        s2ind =  tmp + nz*(j-1 + ny*(i-1));
        
        source_control(s2ind,1) = repmat(l_x(i),nz,1);
        source_control(s2ind,2) = repmat(l_y(j),nz,1);
        source_control(s2ind,3) = l_z(tmp);
        
    end
end

% Plot the source control points
for i = 1:nx * ny * nz
    plot3(source_control(i,1),source_control(i,2),source_control(i,3),'k+');hold on
end

disp('--------- Calculate the source control points mesh (tetrahedrons) -----')
%% Define the mesh for FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global tetra
tetra = [];

for i=1:nx-1
    
    for j=1:ny-1
        
        for k=1:nz-1
            
            ind_tmp =  k + nz*(j-1 + ny*(i-1));
            
            ind_tmp2 =  k+1 + nz*(j-1 + ny*(i-1));
            ind_tmp3 =  k + nz*(j-1 + ny*(i));
            ind_tmp4 =  k + nz*(j + ny*(i));
            ind_tmp5 =  k + nz*(j + ny*(i-1));
            ind_tmp6 =  k+1 + nz*(j + ny*(i));
            ind_tmp7 =  k+1 + nz*(j + ny*(i-1));
            ind_tmp8 =  k+1 + nz*(j-1 + ny*(i));
            
            
            tetra = [tetra;...
                    ind_tmp  ind_tmp2 ind_tmp3 ind_tmp5;...
                    ind_tmp2 ind_tmp3 ind_tmp5 ind_tmp7;...
                    ind_tmp3 ind_tmp7 ind_tmp8 ind_tmp2;...
                    ind_tmp4 ind_tmp7 ind_tmp8 ind_tmp3;...
                    ind_tmp5 ind_tmp7 ind_tmp3 ind_tmp4;...
                    ind_tmp6 ind_tmp7 ind_tmp8 ind_tmp4];
            
            % Just for the plotting
            vari = [source_control(ind_tmp,:);source_control(ind_tmp2,:);source_control(ind_tmp3,:);source_control(ind_tmp4,:);...
                    source_control(ind_tmp5,:);source_control(ind_tmp6,:);source_control(ind_tmp7,:);source_control(ind_tmp8,:)];
            dt = DelaunayTri(vari);
            
            tetramesh(dt);
            
        end
        
    end
    
end

disp('--------- Convert the computed mesh into DelaunayTri class -----')
%% Convert the computed mesh into DelaunayTri class is gonna help us for computing the
%% vertices or tetrahedrons that contain a query of points

global source_tri
global list_edges

trep = TriRep(tetra,source_control);
source_tri = trep;

alpha(.1)
axis equal

list_edges = edges_connected(source_tri);
