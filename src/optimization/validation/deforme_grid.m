%% Generate the deformed MRI

clc
serie = 7;

lava_flex_n      = cellimages2mat(dcmdir.dcmPatient.Study.Series(serie,1).Images);
lava_flex_info = dcmdir.dcmPatient.Study.Series(serie,1).ImagesInfo;

[size1, size2, size3] = size(lava_flex_n);

global lava_flex_ax
global lava_flex_sag
global lava_flex_cor

disp('--------- Image denoising -----')
% Image denoising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size3
    lava_flex(:,:,i) = anisodiff2D(lava_flex_n(:,:,i), 20, 1/7, 30, 1);
end

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
        lava_flex_cor(k,j,1:size1) = lava_flex(1:size1,j,k); % 1:size1 size2-j+1
    end
end

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
    
    tmp = lava_axM{1} * [0 i-1 1]'; % [0 i-1 1], [size1-1 i-1 1]
    [lava_corM{i}, lava_corM_1{i} ,~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 3, tmp(1:3));
    
end


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

disp('--------- Select the slices for each direction -----')

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

disp('--------- Select the slices axial -----')
%% Apply random deformation to the slices
figure;
% Axial
for i = 1:length(slices_ax)
    ind = slices_ax(i);
    
    X_ax_def = [X_ax_def X_ax_v(:,ind)];
    Y_ax_def = [Y_ax_def Y_ax_v(:,ind)];
    Z_ax_def = [Z_ax_def Z_ax_v(:,ind)];
    %% image, M, M1
    vol_ax_eval(:,:,i) = lava_flex_ax(:,:,ind);
    axial_M{i}  = lava_axM{ind}; 
    axial_M1{i} = lava_axM_1{ind}; 

end



% Sagittal
disp('--------- Select the slices sagittal -----')
for i = 1:length(slices_sag)
    ind = slices_sag(i);

    X_sag_def = [X_sag_def X_sag_v(:,ind)];
    Y_sag_def = [Y_sag_def Y_sag_v(:,ind)];
    Z_sag_def = [Z_sag_def Z_sag_v(:,ind)];
    %% image, M, M1
    vol_sag_eval(:,:,i) = lava_flex_sag(:,:,ind);
    sag_M{i}  = lava_sagM{ind}; 
    sag_M1{i} = lava_sagM_1{ind}; 

end


% Coronal
% figure;
disp('--------- Select the slices coronal -----')
for i = 1:length(slices_cor)
    ind = slices_cor(i);

    X_cor_def = [X_cor_def X_cor_v(:,ind)];
    Y_cor_def = [Y_cor_def Y_cor_v(:,ind)];
    Z_cor_def = [Z_cor_def Z_cor_v(:,ind)];
    %% image, M, M1
    vol_cor_eval(:,:,i) = lava_flex_cor(:,:,ind);

    cor_M{i}  = lava_corM{ind}; 
    cor_M1{i} = lava_corM_1{ind}; 

end



%% Calculate the intersection between the different planes
global t 
t = 0:1/14:1;

options = optimset('Display','off');

ax  = slices_ax;
sag = slices_sag;
cor = slices_cor;

global var_cell1_v
global var_array1_v

var_cell1_v  = cell(length(slices_ax),length(slices_sag),length(t));
var_array1_v = zeros(length(slices_ax)*length(slices_sag)*length(t),3);

disp('--------- Calculate the s1 x s2 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(ax)
    
    for j=1:length(sag)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([X_ax_def(:,i)  Y_ax_def(:,i)  Z_ax_def(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([X_sag_def(:,j) Y_sag_def(:,j) Z_sag_def(:,j)]); %% constraints of the sagittal plane
        
        normal = cross([X_sag_def(1,j) Y_sag_def(1,j) Z_sag_def(1,j)]-[X_sag_def(2,j) Y_sag_def(2,j) Z_sag_def(2,j)],[X_sag_def(1,j) Y_sag_def(1,j) Z_sag_def(1,j)]-[X_sag_def(3,j) Y_sag_def(3,j) Z_sag_def(3,j)]); % direction of the intersection line
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        
        [x0,fval,exitflag,output,lambda] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        if exitflag ~= 1 %% means no solution

            for k=1:length(t)
                
                var_cell1_v{i,j,k} = [ ]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell1_v,3) size(var_cell1_v,2) size(var_cell1_v,1)],k,j,i);
                
                var_array1_v(ind_tmp,1) = -Inf;
                var_array1_v(ind_tmp,2) = -Inf;
                var_array1_v(ind_tmp,3) = -Inf;
                
            end
        else
            
            [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
            
            % plot3(V1(:,1),V1(:,2),V1(:,3),'g*');hold on
            
            vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
            
            for k=1:length(t)
                
                plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'r+');hold on;%,'MarkerSize',i);hold on
                var_cell1_v{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell1_v,3) size(var_cell1_v,2) size(var_cell1_v,1)],k,j,i);
                %ind_tmp1 = k + length(t)*(j-1 + length(sag)*(i-1));% it is faster than sub2ind
                
                var_array1_v(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
                var_array1_v(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
                var_array1_v(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
                
            end
        end
        
        
    end
end

disp('--------- Calculate the s1 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global var_cell2_v
global var_array2_v

var_cell2_v = cell(length(slices_ax),length(slices_cor),length(t));
var_array2_v = zeros(length(slices_ax)*length(slices_cor)*length(t),3);

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
                
                var_cell2_v{i,j,k} = []; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell2_v,3) size(var_cell2_v,2) size(var_cell2_v,1)],k,j,i);
                
                var_array2_v(ind_tmp,1) = -Inf;
                var_array2_v(ind_tmp,2) = -Inf;
                var_array2_v(ind_tmp,3) = -Inf;
                
            end
        else
            [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
            
            vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
            
            for k=1:length(t)
                
                plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');hold on%,'MarkerSize',i);hold on
                var_cell2_v{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell2_v,3) size(var_cell2_v,2) size(var_cell2_v,1)],k,j,i);
                %ind_tmp1 = k + length(t)*(j-1 + length(sag)*(i-1));% it is faster than sub2ind
                
                var_array2_v(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
                var_array2_v(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
                var_array2_v(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
                
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

global var_cell3_v
global var_array3_v

var_cell3_v  = cell(length(slices_cor),length(slices_sag),length(t));
var_array3_v = zeros(length(slices_sag)*length(slices_cor)*length(t),3);

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
                
                var_cell3_v{i,j,k} = []; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell3_v,3) size(var_cell3_v,2) size(var_cell3_v,1)],k,j,i);
                
                var_array3_v(ind_tmp,1) = -Inf;
                var_array3_v(ind_tmp,2) = -Inf;
                var_array3_v(ind_tmp,3) = -Inf;
                
            end
        else
            [V1,nr,nre] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!

            vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
            
            for k=1:length(t)
                
                plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'b+'); hold on;
                var_cell3_v{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
                
                ind_tmp = sub2ind([size(var_cell3_v,3) size(var_cell3_v,2) size(var_cell3_v,1)],k,j,i);
                
                var_array3_v(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
                var_array3_v(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
                var_array3_v(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
                
            end
        end
        
        
    end
end

disp('--------- Create the regular grid -----')
%% Create the grid 

n_x = 5;
n_y = 5;
n_z = 3;

% c' {ijk} = c{ijk} + sigma * randn , sigma = deformable amplitude 

X_array = [X_ax_v(:);X_sag_v(:);X_cor_v(:)];
Y_array = [Y_ax_v(:);Y_sag_v(:);Y_cor_v(:)];
Z_array = [Z_ax_v(:);Z_sag_v(:);Z_cor_v(:)];

bb = [min(X_array(:)) max(X_array(:)); ...
      min(Y_array(:)) max(Y_array(:)); ...
      min(Z_array(:)) max(Z_array(:))];

% Create the source control points

l_x = linspace(bb(1,1),bb(1,2),nx);
l_y = linspace(bb(2,1),bb(2,2),ny);
l_z = linspace(bb(3,1),bb(3,2),nz);


control_points = zeros(nx * ny * nz, 3);

for i = 1:nx
    for j = 1:ny
        
        tmp =  1:nz;
        s2ind =  tmp + nz*(j-1 + ny*(i-1));
        
        control_points(s2ind,1) = repmat(l_x(i),nz,1);
        control_points(s2ind,2) = repmat(l_y(j),nz,1);
        control_points(s2ind,3) = l_z(tmp);
        
    end
end

% Plot the source control points
for i = 1:nx * ny * nz
    plot3(control_points(i,1),control_points(i,2),control_points(i,3),'k+');hold on
end

disp('--------- Calculate the source control points mesh (tetra_vhedrons) -----')
%% Define the mesh for FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tetra_c = [];

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
            
            
            tetra_c = [tetra_c;...
                    ind_tmp  ind_tmp2 ind_tmp3 ind_tmp5;...
                    ind_tmp2 ind_tmp3 ind_tmp5 ind_tmp7;...
                    ind_tmp3 ind_tmp7 ind_tmp8 ind_tmp2;...
                    ind_tmp4 ind_tmp7 ind_tmp8 ind_tmp3;...
                    ind_tmp5 ind_tmp7 ind_tmp3 ind_tmp4;...
                    ind_tmp6 ind_tmp7 ind_tmp8 ind_tmp4];
            
            % Just for the plotting
            vari = [control_points(ind_tmp,:);control_points(ind_tmp2,:);control_points(ind_tmp3,:);control_points(ind_tmp4,:);...
                    control_points(ind_tmp5,:);control_points(ind_tmp6,:);control_points(ind_tmp7,:);control_points(ind_tmp8,:)];
            dt = DelaunayTri(vari);
            
            tetramesh(dt);hold on
            
        end
        
    end
    
end
alpha(.1)
axis equal

disp('--------- Convert the computed mesh into DelaunayTri class -----')
%% Convert the computed mesh into DelaunayTri class is gonna help us for computing the
%% vertices or tetra_vhedrons that contain a query of points

trep = TriRep(tetra_c, control_points);
control_tri = trep;

% Deformed grid
sigma = 3; % it is the deformation amplitude
deform_points = control_points + sigma .* randn(size(control_points));

deform_tri = TriRep(tetra_c, deform_points);

%% Plot the two grids
figure;
subplot(121);tetramesh(control_tri);hold on
subplot(122);tetramesh(deform_tri);hold on
            
alpha(.1)
axis equal

disp('--------- Compute the new images -----')
%% Compute the new images


current_tr1 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
current_tr2 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]);
current_tr3 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]);

% New data positions with respect the new computed meshes
b1 = cartToBary(source_tri,  current_tr1,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc1 = baryToCart(deform_tri, current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )

b2 = cartToBary(source_tri,  current_tr2,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc2 = baryToCart(deform_tri, current_tr2, b2); % cartesian coordinates of the points wrt target tri ( p' )

b3 = cartToBary(source_tri,  current_tr3,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc3 = baryToCart(deform_tri, current_tr3, b3); % cartesian coordinates of the points wrt target tri ( p' )


%% Get the intensity values of the p' in the first direction
opt_im_ax  = zeros(size(vol_ax_eval));
opt_im_sag = zeros(size(vol_sag_eval));
opt_im_cor = zeros(size(vol_cor_eval));

r_ax = size(vol_ax_eval,1);
c_ax = size(vol_ax_eval,2);

r_sag = size(vol_sag_eval,1);
c_sag = size(vol_sag_eval,2);

r_cor = size(vol_cor_eval,1);
c_cor = size(vol_cor_eval,2);

for i = 1:size(var_array1_v,1)
    
    [sub_tmp1 sub_tmp2 sub_tmp3] = ind2sub([size(var_cell1_v,3) size(var_cell1_v,2) size(var_cell1_v,1)],i);
    
    tmp_v1 = axial_M1{sub_tmp3} * [xc1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
    
    %% Make sure the indexes are correct and use bilinear interpolation

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r  = min(max(fl2,1),r_ax);
    min_max_r2 = min(max(cl2,1),r_ax);
    min_max_c  = min(max(fl,1), c_ax);
    min_max_c2 = min(max(cl,1), c_ax);
    
    neig = [vol_ax_eval(min_max_r, min_max_c, sub_tmp3(i))   vol_ax_eval(min_max_r, min_max_c2, sub_tmp3(i));...
            vol_ax_eval(min_max_r2,min_max_c, sub_tmp3(i))   vol_ax_eval(min_max_r2,min_max_c2, sub_tmp3(i))];
            
    new_i = min_max_r;
    new_j = min_max_c;
    new_k = sub_tmp3;

    opt_im_ax(new_i, new_j, new_k) = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));


end 

for i = 1:size(var_array1_v,1)
    
    [sub_tmp1 sub_tmp2 sub_tmp3] = ind2sub([size(var_cell1_v,3) size(var_cell1_v,2) size(var_cell1_v,1)],i);
    
    tmp_v1 = sag_M1{sub_tmp2} * [xc1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
    
    %% Make sure the indexes are correct and use bilinear interpolation

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r  = min(max(fl2,1),r_sag);
    min_max_r2 = min(max(cl2,1),r_sag);
    min_max_c  = min(max(fl,1), c_sag);
    min_max_c2 = min(max(cl,1), c_sag);
    
    neig = [vol_sag_eval(min_max_r, min_max_c, sub_tmp2(i))   vol_sag_eval(min_max_r, min_max_c2, sub_tmp2(i));...
            vol_sag_eval(min_max_r2,min_max_c, sub_tmp2(i))   vol_sag_eval(min_max_r2,min_max_c2, sub_tmp2(i))];
            
    new_i = min_max_r;
    new_j = min_max_c;
    new_k = sub_tmp2;

    opt_im_sag(new_i, new_j, new_k) = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));


end 

for i = 1:size(var_array3_v,1)
    
    [sub_tmp1 sub_tmp2 sub_tmp3] = ind2sub([size(var_cell3_v,3) size(var_cell3_v,2) size(var_cell3_v,1)],i);
    
    tmp_v1 = cor_M1{sub_tmp3} * [xc3(i,:) 1]'; % 3D point to 2D point in the frame coordinates
    
    %% Make sure the indexes are correct and use bilinear interpolation

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r  = min(max(fl2,1),r_sag);
    min_max_r2 = min(max(cl2,1),r_sag);
    min_max_c  = min(max(fl,1), c_sag);
    min_max_c2 = min(max(cl,1), c_sag);
    
    neig = [vol_sag_eval(min_max_r, min_max_c, sub_tmp3(i))   vol_sag_eval(min_max_r, min_max_c2, sub_tmp3(i));...
            vol_sag_eval(min_max_r2,min_max_c, sub_tmp3(i))   vol_sag_eval(min_max_r2,min_max_c2, sub_tmp3(i))];
            
    new_i = min_max_r;
    new_j = min_max_c;
    new_k = sub_tmp3;

    opt_im_cor(new_i, new_j, new_k) = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));


end 
    
figure;subplot(121);imshow( vol_ax_eval(:,:,9),[]); subplot(122);imshow( opt_im_ax(:,:,9),[]);
figure;subplot(121);imshow(vol_sag_eval(:,:,10),[]); subplot(122);imshow(opt_im_sag(:,:,10),[]);
figure;subplot(121);imshow(vol_cor_eval(:,:,10),[]); subplot(122);imshow(opt_im_cor(:,:,10),[]);


