%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read the LAVA-Flex in the 3 directions
%% and apply a random rigid transformation
%% to each one
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
serie = 7;

lava_flex      = cellimages2mat(dcmdir.dcmPatient.Study.Series(serie,1).Images);
lava_flex_info = dcmdir.dcmPatient.Study.Series(serie,1).ImagesInfo;

[size1, size2, size3] = size(lava_flex);

global lava_flex_ax
global lava_flex_sag
global lava_flex_cor

disp('--------- Image denoising -----')
% Image denoising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:size3
%     lava_flex(:,:,i) = anisodiff2D(lava_flex_n(:,:,i), 20, 1/7, 30, 1);
% end

disp('--------- Define each direction -----')
lava_flex_ax  = lava_flex;

for k = 1:size3
    for i = 1:size1
        lava_flex_sag(k,i,:) = lava_flex(i,:,k); %size1-i+1,size2:-1:1,k
%         lava_flex_sag(k,i,1:size2) = lava_flex(i,1:size2,k); %size1-i+1,size2:-1:1,k
    end
end


for k = 1:size3
    for j = 1:size2
        lava_flex_cor(k,j,1:size1) = lava_flex(1:size1,size2-j+1,k); % 1:size1 size2-j+1
    end
end

% for k = 1:size3
%     for i = 1:size1
%         lava_flex_sag(k,i,1:size2) = lava_flex(i,size2:-1:1,k); %size1-i+1,size2:-1:1,k
% %         lava_flex_sag(k,i,1:size2) = lava_flex(i,1:size2,k); %size1-i+1,size2:-1:1,k
%     end
% end

%% Apply random transform to each view, in each plane
%% First try with rigid transforms
global lava_axM
global lava_axM_1
global lava_sagM
global lava_sagM_1
global lava_corM
global lava_corM_1

% Translation
global ortho
ortho = 1;

disp('--------- Calculate transformations -----')
%% First compute the transformation M, from 2D image to the 3D RCS
for i = 1:size3
    
    [lava_axM{i},lava_axM_1{i},~] = compute_M_M1(lava_flex_info{i}, ortho, 1);
    
end

for i = 1:size2
    
    tmp = lava_axM{1} * [i-1 0 1]'; % size2-i+1, [i-1 0 1] [size2-i 0 1]
    [lava_sagM{i}, lava_sagM_1{i} ,~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 2, tmp(1:3));
    
end

for i = 1:size1
    
    tmp = lava_axM{1} * [size1-1 i-1 1]'; % [0 i-1 1], [size1-1 i-1 1]
    [lava_corM{i}, lava_corM_1{i} ,~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 3, tmp(1:3));
    
end

lava_axM{1} * [size2-1 0 1]'
lava_axM{1} * [size1-1 0 1]'

X_ax_v = [];
Y_ax_v = [];
Z_ax_v = [];

X_sag_v = [];
Y_sag_v = [];
Z_sag_v = [];

X_cor_v = [];
Y_cor_v = [];
Z_cor_v = [];


%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('Display','off');

disp('--------- Calculate M^(-1) and plane eq. for each slice in the first direction -----')

figure;
for ax = 1:size(lava_flex_ax,3)
    
    [x,y,z] = calculate4corners( lava_axM{ax} );
    
    X_ax_v = [X_ax_v x'];
    Y_ax_v = [Y_ax_v y'];
    Z_ax_v = [Z_ax_v z'];
    
    points = [x' y' z'];
    plot3(points(:,1),points(:,2),points(:,3),'b+');hold on
    
    if ax == 1
        N1 = cross([X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(2,ax) Y_ax_v(2,ax) Z_ax_v(2,ax)],[X_ax_v(1,ax) Y_ax_v(1,ax) Z_ax_v(1,ax)]-[X_ax_v(3,ax) Y_ax_v(3,ax) Z_ax_v(3,ax)]); % normal to the axial (ax)
        N1 = N1./norm(N1)
    end
    
end

disp('--------- Calculate M^(-1) and plane eq. for each slice in the third direction -----')

for cor = 1:size1
    
%     [x,y,z] = calculate4corners_cor( lava_axM{1}, lava_axM{end}, [0 cor-1;511 cor-1]); %  [cor-1 511] , [cor-1 0] [0 cor-1;511 cor-1]
    [x,y,z] = calculate4corners( lava_corM{cor} );
    
    X_cor_v = [X_cor_v x'];
    Y_cor_v = [Y_cor_v y'];
    Z_cor_v = [Z_cor_v z'];
    
    points = [x' y' z'];
    plot3(points(:,1),points(:,2),points(:,3),'r*');hold on
    
    if cor == 1
        N3 = cross([X_cor_v(1,cor) Y_cor_v(1,cor) Z_cor_v(1,cor)]-[X_cor_v(2,cor) Y_cor_v(2,cor) Z_cor_v(2,cor)],[X_cor_v(1,cor) Y_cor_v(1,cor) Z_cor_v(1,cor)]-[X_cor_v(3,cor) Y_cor_v(3,cor) Z_cor_v(3,cor)]); % normal to the axial (ax)
        N3 = N3./norm(N3)
    end
    
end

disp('--------- Calculate M^(-1) and plane eq. for each slice in the second direction -----')

for sag = 1:size2
    
    %     [x,y,z] = calculate4corners_cor( lava_axM{1}, lava_axM{end}, [511-sag 511;511-sag 0]); %  [0 511-sag;511 511-sag]
    [x,y,z] = calculate4corners( lava_sagM{sag} );
    
    X_sag_v = [X_sag_v x'];
    Y_sag_v = [Y_sag_v y'];
    Z_sag_v = [Z_sag_v z'];
    
    points = [x' y' z'];
    plot3(points(:,1),points(:,2),points(:,3),'g*');hold on
    
    if sag == 1
        N2 = cross([X_sag_v(1,sag) Y_sag_v(1,sag) Z_sag_v(1,sag)]-[X_sag_v(2,sag) Y_sag_v(2,sag) Z_sag_v(2,sag)],[X_sag_v(1,sag) Y_sag_v(1,sag) Z_sag_v(1,sag)]-[X_sag_v(3,sag) Y_sag_v(3,sag) Z_sag_v(3,sag)]); % normal to the axial (ax)
        N2 = N2./norm(N2)
    end
    
end

%% Select a few slices, like in T2W scans
n_slices_ax = 25;
slices_ax = 4:round(size3/n_slices_ax):size3-4;

st_sag   = 100;
end_sag = size2 - st_sag;
n_slices_sag = 25;
slices_sag = st_sag:ceil((end_sag-st_sag)/n_slices_sag):end_sag;

st_cor   = 150;
end_cor = size1 - st_cor;
n_slices_cor = 25;
slices_cor = st_cor:ceil((end_cor-st_cor)/n_slices_cor):end_cor;

X_ax_def = [];
Y_ax_def = [];
Z_ax_def = [];

X_sag_def = [];
Y_sag_def = [];
Z_sag_def = [];

X_cor_def = [];
Y_cor_def = [];
Z_cor_def = [];

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global axial_M
global axial_M1
global sag_M
global sag_M1
global cor_M
global cor_M1

%% Apply random deformation to the slices
% Axial
for i = 1:length(slices_ax)
    ind = slices_ax(i);
    [out_plane, tr] = apply_t([X_ax_v(:,ind) Y_ax_v(:,ind) Z_ax_v(:,ind)],[7 7 0]');
    X_ax_def = [X_ax_def out_plane(:,1)];
    Y_ax_def = [Y_ax_def out_plane(:,2)];
    Z_ax_def = [Z_ax_def out_plane(:,3)];
    %% image, M, M1
    vol_ax_eval(:,:,i) = lava_flex_ax(:,:,ind);
    axial_M{i}  = tr * lava_axM{ind};
    axial_M1{i} = lava_axM_1{ind} / tr;
    %     plot3(X_ax_v(1,ind), Y_ax_v(1,ind), Z_ax_v(1,ind),'m+');hold on
    %     plot3(X_ax_def(1,1), Y_ax_def(1,1), Z_ax_def(1,1), 'g*');hold on
end



% Sagittal

for i = 1:length(slices_sag)
    ind = slices_sag(i);
    [out_plane, tr] = apply_t([X_sag_v(:,ind) Y_sag_v(:,ind) Z_sag_v(:,ind)],[0 7 2]');
    X_sag_def = [X_sag_def out_plane(:,1)];
    Y_sag_def = [Y_sag_def out_plane(:,2)];
    Z_sag_def = [Z_sag_def out_plane(:,3)];
    %% image, M, M1
    vol_sag_eval(:,:,i) = lava_flex_sag(:,:,ind);
    sag_M{i}  = tr * lava_sagM{ind};
    sag_M1{i} = lava_sagM_1{ind} / tr;
    %     plot3(X_sag_v(1,ind), Y_sag_v(1,ind), Z_sag_v(1,ind),'m+');hold on
    %     plot3(X_sag_def(1,1), Y_sag_def(1,1), Z_sag_def(1,1), 'g*');hold on
end


% Coronal
% figure;
for i = 1:length(slices_cor)
    ind = slices_cor(i);
    [out_plane, tr] = apply_t([X_cor_v(:,ind) Y_cor_v(:,ind) Z_cor_v(:,ind)],[-5 0 7]');
    X_cor_def = [X_cor_def out_plane(:,1)];
    Y_cor_def = [Y_cor_def out_plane(:,2)];
    Z_cor_def = [Z_cor_def out_plane(:,3)];
    %% image, M, M1
    vol_cor_eval(:,:,i) = lava_flex_cor(:,:,ind);
    cor_M{i}  = tr * lava_corM{ind};
    cor_M1{i} = lava_corM_1{ind} / tr;
    %     plot3(X_sag_v(1,ind), Y_sag_v(1,ind), Z_sag_v(1,ind),'m+');hold on
    %     plot3(X_sag_def(1,1), Y_sag_def(1,1), Z_sag_def(1,1), 'g*');hold on
end

% figure;
% fill3(X_ax_def(:,1),Y_ax_def(:,1),Z_ax_def(:,1),'r');hold on % first axial plane
% fill3(X_ax_def(:,length(slices_ax)),Y_ax_def(:,length(slices_ax)),Z_ax_def(:,length(slices_ax)),'r');hold on % first axial plane
% fill3(X_sag_def(:,1),Y_sag_def(:,1),Z_sag_def(:,1),'b');hold on % first sagittal plane
% fill3(X_sag_def(:,length(slices_sag)),Y_sag_def(:,length(slices_sag)),Z_sag_def(:,length(slices_sag)),'b');hold on % second sagittal plane
% alpha(.2)

%% Calculate the intersection between the different planes
global t 
t = 0:1/14:1;

options = optimset('Display','off');

ax  = slices_ax;
sag = slices_sag;
cor = slices_cor;

global var_cell1
global var_array1

var_cell1  = cell(length(slices_ax),length(slices_sag),length(t));
var_array1 = zeros(length(slices_ax)*length(slices_sag)*length(t),3);

disp('--------- Calculate the s1 x s2 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(ax)
    
    for j=1:length(sag)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([X_ax_def(:,i) Y_ax_def(:,i) Z_ax_def(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([X_sag_def(:,j) Y_sag_def(:,j) Z_sag_def(:,j)]); %% constraints of the sagittal plane
        
        normal = cross([X_sag_def(1,j) Y_sag_def(1,j) Z_sag_def(1,j)]-[X_sag_def(2,j) Y_sag_def(2,j) Z_sag_def(2,j)],[X_sag_def(1,j) Y_sag_def(1,j) Z_sag_def(1,j)]-[X_sag_def(3,j) Y_sag_def(3,j) Z_sag_def(3,j)]); % direction of the intersection line
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        
        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        if exitflag ~= 1 %% means no solution
            'hola1'
            for k=1:length(t)
                
                var_cell1{i,j,k} = [ ]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell1,3) size(var_cell1,2) size(var_cell1,1)],k,j,i);
                
                var_array1(ind_tmp,1) = -Inf;
                var_array1(ind_tmp,2) = -Inf;
                var_array1(ind_tmp,3) = -Inf;
                
            end
        else
            
            [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
            
            % plot3(V1(:,1),V1(:,2),V1(:,3),'g*');hold on
            
            vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
            
            for k=1:length(t)
                
                plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'r+');hold on;%,'MarkerSize',i);hold on
                var_cell1{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell1,3) size(var_cell1,2) size(var_cell1,1)],k,j,i);
                %ind_tmp1 = k + length(t)*(j-1 + length(sag)*(i-1));% it is faster than sub2ind
                
                var_array1(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
                var_array1(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
                var_array1(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
                
            end
        end
        
        
    end
end

disp('--------- Calculate the s1 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global var_cell2
global var_array2

var_cell2 = cell(length(slices_ax),length(slices_cor),length(t));
var_array2 = zeros(length(slices_ax)*length(slices_cor)*length(t),3);

options = optimset('Display','off');


for i=1:length(ax)
    
    for j=1:length(cor)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([X_ax_def(:,i) Y_ax_def(:,i) Z_ax_def(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([X_cor_def(:,j) Y_cor_def(:,j) Z_cor_def(:,j)]); %% constraints of the sagittal plane
        
        if i==1 && j==1
            normal = cross([X_cor_def(1,j) Y_cor_def(1,j) Z_cor_def(1,j)]-[X_cor_def(2,j) Y_cor_def(2,j) Z_cor_def(2,j)],[X_cor_def(1,j) Y_cor_def(1,j) Z_cor_def(1,j)]-[X_cor_def(3,j) Y_cor_def(3,j) Z_cor_def(3,j)]); % direction of the intersection line
        end
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        
        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        if exitflag ~= 1 %% means no solution
            'hola2'
            for k=1:length(t)
                
                var_cell2{i,j,k} = []; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell2,3) size(var_cell2,2) size(var_cell2,1)],k,j,i);
                
                var_array2(ind_tmp,1) = -Inf;
                var_array2(ind_tmp,2) = -Inf;
                var_array2(ind_tmp,3) = -Inf;
                
            end
        else
            [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
            
            vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
            
            for k=1:length(t)
                
                plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');hold on%,'MarkerSize',i);hold on
                var_cell2{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell2,3) size(var_cell2,2) size(var_cell2,1)],k,j,i);
                %ind_tmp1 = k + length(t)*(j-1 + length(sag)*(i-1));% it is faster than sub2ind
                
                var_array2(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
                var_array2(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
                var_array2(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
                
            end
        end
        
        
    end
end


% fill3(X_ax(:,1),Y_ax(:,1),Z_ax(:,1),'r');hold on % first axial plane
% fill3(X_ax(:,total_ax),Y_ax(:,total_ax),Z_ax(:,total_ax),'r');hold on % first axial plane
% fill3(X_cor(:,1),Y_cor(:,1),Z_cor(:,1),'b');hold on % first sagittal plane
% fill3(X_cor(:,total_cor),Y_cor(:,total_cor),Z_cor(:,total_cor),'b');hold on % second sagittal plane
% alpha(.2)
disp('--------- Calculate the s2 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global var_cell3
global var_array3

var_cell3  = cell(length(slices_cor),length(slices_sag),length(t));
var_array3 = zeros(length(slices_sag)*length(slices_cor)*length(t),3);

options = optimset('Display','off');

for i=1:length(cor)
    
    for j=1:length(sag)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([X_cor_def(:,i) Y_cor_def(:,i) Z_cor_def(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([X_sag_def(:,j) Y_sag_def(:,j) Z_sag_def(:,j)]); %% constraints of the sagittal plane
        
        if i==1 && j==1
            normal = cross([X_sag_def(1,j) Y_sag_def(1,j) Z_sag_def(1,j)]-[X_sag_def(2,j) Y_sag_def(2,j) Z_sag_def(2,j)],[X_sag_def(1,j) Y_sag_def(1,j) Z_sag_def(1,j)]-[X_sag_def(3,j) Y_sag_def(3,j) Z_sag_def(3,j)]); % direction of the intersection line
        end
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        
        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        if exitflag  ~= 1
            'hola3'
            for k=1:length(t)
                
                var_cell3{i,j,k} = []; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell3,3) size(var_cell3,2) size(var_cell3,1)],k,j,i);
                
                var_array3(ind_tmp,1) = -Inf;
                var_array3(ind_tmp,2) = -Inf;
                var_array3(ind_tmp,3) = -Inf;
                
            end
        else
            [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!

            vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
            
            for k=1:length(t)
                
                plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'b+'); hold on;
                var_cell3{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell3,3) size(var_cell3,2) size(var_cell3,1)],k,j,i);
                
                var_array3(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
                var_array3(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
                var_array3(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
                
            end
        end
        
        
    end
end

disp('--------- Calculate the source control points ( # N^3 ) -----')
%% Calculate the bounding box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bb = [xmin xmax;ymin ymax;zmin zmax]
global var_array

var_array = [var_array1;var_array2;var_array3];

tmp_var1 = var_array(:,1);
tmp_var2 = var_array(:,2);
tmp_var3 = var_array(:,3);

bb = [min(min(tmp_var1(tmp_var1~=-Inf))) max(max(var_array(:,1))); ...
      min(min(tmp_var2(tmp_var2~=-Inf))) max(max(var_array(:,2))); ...
      min(min(tmp_var3(tmp_var3~=-Inf))) max(max(var_array(:,3)))];

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
figure;
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
            
            tetramesh(dt);hold on
            
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

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------- Start Optimization  -----')

disp('--------- Preparing first intersection: axial/sagittal  -----')
size1_int1 = size(var_cell1,3);
size2_int1 = size(var_cell1,2);
size3_int1 = size(var_cell1,1);

global sub_1_int1
global sub_2_int1
global sub_3_int1
global current_tr_int1
global c2b_coord_int1


[sub_1_int1, sub_2_int1, sub_3_int1] = ind2sub([size1_int1 size2_int1 size3_int1],1:size(var_array1,1));
current_tr_int1 = tsearchn(source_tri.X,source_tri.Triangulation,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int1  = cartToBary(source_tri,current_tr_int1,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % barycentric coordinates of the points wrt source tri

disp('--------- Preparing second intersection: axial/coronal  -----')
size1_int2 = size(var_cell2,3);
size2_int2 = size(var_cell2,2);
size3_int2 = size(var_cell2,1);

global sub_1_int2
global sub_2_int2
global sub_3_int2
global current_tr_int2
global c2b_coord_int2


[sub_1_int2, sub_2_int2, sub_3_int2] = ind2sub([size1_int2 size2_int2 size3_int2],1:size(var_array2,1));
current_tr_int2 = tsearchn(source_tri.X,source_tri.Triangulation,[var_array2(:,1) var_array2(:,2) var_array2(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int2  = cartToBary(source_tri,current_tr_int2,[var_array2(:,1) var_array2(:,2) var_array2(:,3)]); % barycentric coordinates of the points wrt source tri


disp('--------- Preparing third intersection: sagittal/coronal  -----')
size1_int3 = size(var_cell3,3);
size2_int3 = size(var_cell3,2);
size3_int3 = size(var_cell3,1);

global sub_1_int3
global sub_2_int3
global sub_3_int3
global current_tr_int3
global c2b_coord_int3


[sub_1_int3, sub_2_int3, sub_3_int3] = ind2sub([size1_int3 size2_int3 size3_int3],1:size(var_array3,1));
current_tr_int3 = tsearchn(source_tri.X,source_tri.Triangulation,[var_array3(:,1) var_array3(:,2) var_array3(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int3  = cartToBary(source_tri,current_tr_int3,[var_array3(:,1) var_array3(:,2) var_array3(:,3)]); % barycentric coordinates of the points wrt source tri


lbv = [-20 -20 -20];
ubv = [ 20  20  20];

% lb = repmat(lbv,2*size(source_tri.X,1),1) + [source_tri.X;source_tri.X];
% ub = repmat(ubv,2*size(source_tri.X,1),1) + [source_tri.X;source_tri.X];
lb = [];
ub = [];

opts = optimset('Jacobian','on','Display','iter','MaxIter', 20);

tic

preparing_eval;

tim = toc

tic
%% initialization
mesh0 = [source_tri.X(:,1:2);source_tri.X(:,2:3);source_tri.X(:,1) source_tri.X(:,3)]; % [source_tri.X(:,1:2);source_tri.X(:,2:3)]

options = optimset('Display','iter','MaxIter', 40);
[xfinal_tmp fval exitflag output] = fminunc(@myfun_unc_ortho_eval, mesh0, options);

xfinal(1:size(source_tri.X,1),:) = [xfinal_tmp(1:size(source_tri.X,1),:) source_tri.X(:,3)];
xfinal(1+size(source_tri.X,1):2*size(source_tri.X,1),:)   = [source_tri.X(:,1) xfinal_tmp(1+size(source_tri.X,1):2*size(source_tri.X,1),:)];
xfinal(1+2*size(source_tri.X,1):3*size(source_tri.X,1),:) = [ xfinal_tmp(1+2*size(source_tri.X,1):3*size(source_tri.X,1),1) source_tri.X(:,2) xfinal_tmp(1+2*size(source_tri.X,1):3*size(source_tri.X,1),2)];

%% Define the initial mesh as a vector
% mesh_vector_tmp = [source_tri.X;source_tri.X]';
% mesh_vector = mesh_vector_tmp(:);
% 
% [xfinal_v, fval] = optimizationWithDE(1, 6*size(source_tri.X,1),[],[],[], [], [], [], [], []);
% xfinal = reshape(xfinal_v',2*size(source_tri.X,1),3) +  [source_tri.X;source_tri.X];

tim = toc