%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read the LAVA-Flex in the 3 directions
%% and apply a random rigid transformation
%% to each one
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
serie = 7;

save_name = 'deformation_0.mat';

lava_flex_n      = cellimages2mat(dcmdir.dcmPatient.Study.Series(serie,1).Images);
lava_flex_info = dcmdir.dcmPatient.Study.Series(serie,1).ImagesInfo;

[size1, size2, size3] = size(lava_flex_n);

global lava_flex_ax
global lava_flex_sag
global lava_flex_cor

disp('--------- Image denoising -----')
% lava_flex = lava_flex_n;
% Image denoising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size3
    lava_flex(:,:,i) = anisodiff2D(lava_flex_n(:,:,i), 20, 1/7, 30, 1);
end

disp('--------- Define each direction -----')
lava_flex_ax  = lava_flex;

for k = 1:size3
    for i = 1:size1
        lava_flex_sag(k,i,:) = lava_flex(i,:,k); 
    end
end


for k = 1:size3
    for j = 1:size2
        lava_flex_cor(k,j,:) = lava_flex(:,j,k); 
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
    
    [lava_axM{i},lava_axM_1{i}, ~] = compute_M_M1(lava_flex_info{i}, ortho, 1);
    
end


for i = 1:size2
    
    tmp = lava_axM{1} * [i-1 0 1]';
    [lava_sagM{i}, lava_sagM_1{i} , ~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 2, tmp(1:3));
    
end

for i = 1:size1
    
    tmp = lava_axM{1} * [0 i-1 1]';
    [lava_corM{i}, lava_corM_1{i} , ~] = compute_M_M1_synt(lava_flex_info{1}, ortho, 3, tmp(1:3));
    
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
    
    [x,y,z] = calculate4corners( lava_axM{ax}, [0 size2-1], [0 size1-1] );
    
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
    [x,y,z] = calculate4corners( lava_corM{cor}, [0 size2-1], [0 size3-1] );
    
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
    [x,y,z] = calculate4corners( lava_sagM{sag}, [0 size1-1], [0 size3-1]  );
    
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

disp('--------- Apply deformations -----')

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

% X_ax_def = X_ax_v;
% Y_ax_def = Y_ax_v;
% Z_ax_def = Z_ax_v;
% 
% X_sag_def = X_sag_v;
% Y_sag_def = Y_sag_v;
% Z_sag_def = Z_sag_v;
% 
% X_cor_def = X_cor_v;
% Y_cor_def = Y_cor_v;
% Z_cor_def = Z_cor_v;

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global axial_M
global axial_M1
global sag_M
global sag_M1
global cor_M
global cor_M1

disp('--------- Apply deformations axial -----')
%% Apply random deformation to the slices
figure;
% Axial
for i = 1:length(slices_ax)
    ind = slices_ax(i);
%     [out_plane, tr] = apply_t([X_ax_v(:,ind) Y_ax_v(:,ind) Z_ax_v(:,ind)],[7 7 0]');
%     [T1, VX_ax(:,:,i), VY_ax(:,:,i)] = solving_ad3( double(lava_flex_ax(:,:,ind)) );
%      X_ax_def = [X_ax_def out_plane(:,1)];
%      Y_ax_def = [Y_ax_def out_plane(:,2)];
%      Z_ax_def = [Z_ax_def out_plane(:,3)];
    X_ax_def = [X_ax_def X_ax_v(:,ind)];
    Y_ax_def = [Y_ax_def Y_ax_v(:,ind)];
    Z_ax_def = [Z_ax_def Z_ax_v(:,ind)];
    %% image, M, M1
    vol_ax_eval(:,:,i) = lava_flex_ax(:,:,ind); % opt_im_ax(:,:,i);%
%     vol_ax_eval(:,:,i) = T1;
    axial_M{i}  = lava_axM{ind}; % tr * 
    axial_M1{i} = lava_axM_1{ind}; %  / tr
%         plot3(X_ax_v(1,ind), Y_ax_v(1,ind), Z_ax_v(1,ind),'m+');hold on
%         plot3(X_ax_def(1,1), Y_ax_def(1,1), Z_ax_def(1,1), 'g*');hold on
end



% Sagittal
disp('--------- Apply deformations sagittal -----')
for i = 1:length(slices_sag)
    ind = slices_sag(i);
%     [out_plane, tr] = apply_t([X_sag_v(:,ind) Y_sag_v(:,ind) Z_sag_v(:,ind)],[0 0 0]');
%     [T1, VX_sag(:,:,i), VY_sag(:,:,i)] = solving_ad3( double(lava_flex_sag(:,:,ind)) );
%     X_sag_def = [X_sag_def out_plane(:,1)];
%     Y_sag_def = [Y_sag_def out_plane(:,2)];
%     Z_sag_def = [Z_sag_def out_plane(:,3)];
    X_sag_def = [X_sag_def X_sag_v(:,ind)];
    Y_sag_def = [Y_sag_def Y_sag_v(:,ind)];
    Z_sag_def = [Z_sag_def Z_sag_v(:,ind)];
    %% image, M, M1
    vol_sag_eval(:,:,i) = lava_flex_sag(:,:,ind); %opt_im_sag(:,:,i);%
%     vol_sag_eval(:,:,i) = T1;
    sag_M{i}  = lava_sagM{ind}; % tr * 
    sag_M1{i} = lava_sagM_1{ind}; %  / tr
%         plot3(X_sag_v(1,ind), Y_sag_v(1,ind), Z_sag_v(1,ind),'m+');hold on
%         plot3(X_sag_def(1,1), Y_sag_def(1,1), Z_sag_def(1,1), 'g*');hold on
end


% Coronal
% figure;
disp('--------- Apply deformations coronal -----')
for i = 1:length(slices_cor)
    ind = slices_cor(i);
%     [out_plane, tr] = apply_t([X_cor_v(:,ind) Y_cor_v(:,ind) Z_cor_v(:,ind)],[0 0 0]');
%     [T1, VX_cor(:,:,i), VY_cor(:,:,i)] = solving_ad3( double(lava_flex_cor(:,:,ind)) );
    
%     X_cor_def = [X_cor_def out_plane(:,1)];
%     Y_cor_def = [Y_cor_def out_plane(:,2)];
%     Z_cor_def = [Z_cor_def out_plane(:,3)];
    X_cor_def = [X_cor_def X_cor_v(:,ind)];
    Y_cor_def = [Y_cor_def Y_cor_v(:,ind)];
    Z_cor_def = [Z_cor_def Z_cor_v(:,ind)];
    %% image, M, M1
    vol_cor_eval(:,:,i) = lava_flex_cor(:,:,ind); %opt_im_cor(:,:,i);%
%     vol_cor_eval(:,:,i) = T1;
    cor_M{i}  =  lava_corM{ind}; % tr * 
    cor_M1{i} = lava_corM_1{ind}; %  / tr
%         plot3(X_sag_v(1,ind), Y_sag_v(1,ind), Z_sag_v(1,ind),'m+');hold on
%         plot3(X_sag_def(1,1), Y_sag_def(1,1), Z_sag_def(1,1), 'g*');hold on
end

% figure;
% fill3(X_ax_def(:,1),Y_ax_def(:,1),Z_ax_def(:,1),'r');hold on % first axial plane
% fill3(X_ax_def(:,length(slices_ax)),Y_ax_def(:,length(slices_ax)),Z_ax_def(:,length(slices_ax)),'r');hold on % first axial plane
% fill3(X_sag_def(:,1),Y_sag_def(:,1),Z_sag_def(:,1),'b');hold on % first sagittal plane
% fill3(X_sag_def(:,length(slices_sag)),Y_sag_def(:,length(slices_sag)),Z_sag_def(:,length(slices_sag)),'b');hold on % second sagittal plane
% alpha(.2)

%% Calculate the intersection between the different planes
global t 
t = 0:1/34:1;

options = optimset('Display','off');

ax  = slices_ax;
sag = slices_sag;
cor = slices_cor;

global var_cell1_v
global var_array1_v

var_cell1_v  = cell(length(slices_ax),length(slices_sag),length(t));
var_array1_v = zeros(length(slices_ax)*length(slices_sag)*length(t),3);

new_im_ax = zeros(size(vol_ax_eval));
new_im_sag = zeros(size(vol_sag_eval));

new_im_ax2 = zeros(size(vol_ax_eval));
new_im_cor1 = zeros(size(vol_cor_eval));

new_im_cor2 = zeros(size(vol_cor_eval));
new_im_sag2 = zeros(size(vol_sag_eval));

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
            'hola1'
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
                
                %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [i1, j1, i2, j2, real_v]  = compute_coord(axial_M1{i}, [var_cell1_v{i,j,k} 1], size(vol_ax_eval,1), size(vol_ax_eval,2));
                
                neig = [vol_ax_eval(i1, j1, i)   vol_ax_eval(i1, j2, i);...
                        vol_ax_eval(i2, j1, i)   vol_ax_eval(i2, j2, i)];
                
                new_im_ax(i1, j1, i) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
                diff1 = new_im_ax(i1, j1, i);
                %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [i1, j1, i2, j2, real_v]  = compute_coord(sag_M1{j}, [var_cell1_v{i,j,k} 1], size(vol_sag_eval,1), size(vol_sag_eval,2));
                
                neig = [vol_sag_eval(i1, j1, j)   vol_sag_eval(i1, j2, j);...
                        vol_sag_eval(i2, j1, j)   vol_sag_eval(i2, j2, j)];
                
                new_im_sag(i1, j1, j) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
                diff2 = new_im_sag(i1, j1, j);
                
                diff(ind_tmp) = abs(diff1 - diff2);
            
            end
        end
        
        
    end
end
% show_results(new_im_ax);
% % show_results(new_im_sag);

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
                
                %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [i1, j1, i2, j2, real_v]  = compute_coord(axial_M1{i}, [var_cell2_v{i,j,k} 1], size(vol_ax_eval,1), size(vol_ax_eval,2));
                
                neig = [vol_ax_eval(i1, j1, i)   vol_ax_eval(i1, j2, i);...
                        vol_ax_eval(i2, j1, i)   vol_ax_eval(i2, j2, i)];
                
                new_im_ax2(i1, j1, i) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
                diff1 = new_im_ax(i1, j1, i);
                %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [i1, j1, i2, j2, real_v]  = compute_coord(cor_M1{j}, [var_cell2_v{i,j,k} 1], size(vol_cor_eval,1), size(vol_cor_eval,2));

                neig = [vol_cor_eval(i1, j1, j)   vol_cor_eval(i1, j2, j);...
                        vol_cor_eval(i2, j1, j)   vol_cor_eval(i2, j2, j)];
                
                new_im_cor1(i1, j1, j) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
                diff2 = new_im_cor1(i1, j1, j);
                
                diff(ind_tmp) = abs(diff1 - diff2);
                
            end
        end
        
        
    end
end

% show_results(new_im_ax2);
% show_results(new_im_cor1);

% figure;
% fill3(X_ax_def(:,1),Y_ax_def(:,1),Z_ax_def(:,1),'r');hold on % first axial plane
fill3(X_ax_def(:,i),Y_ax_def(:,i),Z_ax_def(:,i),'r');hold on % first axial plane
% fill3(X_cor_def(:,1),Y_cor_def(:,1),Z_cor_def(:,1),'b');hold on % first sagittal plane
fill3(X_cor_def(:,j),Y_cor_def(:,j),Z_cor_def(:,j),'b');hold on % second sagittal plane
alpha(.2)
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
                
                %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [i1, j1, i2, j2, real_v]  = compute_coord(cor_M1{i}, [var_cell3_v{i,j,k} 1], size(vol_cor_eval,1), size(vol_cor_eval,2));
                
                neig = [vol_cor_eval(i1, j1, i)   vol_cor_eval(i1, j2, i);...
                        vol_cor_eval(i2, j1, i)   vol_cor_eval(i2, j2, i)];
                
                new_im_cor2(i1, j1, i) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
                diff1 = new_im_cor2(i1, j1, i);
                %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [i1, j1, i2, j2, real_v]  = compute_coord(sag_M1{j}, [var_cell3_v{i,j,k} 1], size(vol_sag_eval,1), size(vol_sag_eval,2));
                
                neig = [vol_sag_eval(i1, j1, j)   vol_sag_eval(i1, j2, j);...
                        vol_sag_eval(i2, j1, j)   vol_sag_eval(i2, j2, j)];
                
                new_im_sag2(i1, j1, j) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
                diff2 = new_im_sag2(i1, j1, j);
                
                diff(ind_tmp) = abs(diff1 - diff2);
                
            end
        end
        
        
    end
end
show_results(new_im_cor2);
show_results(new_im_sag2);

% figure;
% fill3(X_sag_def(:,1),Y_sag_def(:,1),Z_sag_def(:,1),'r');hold on % first axial plane
fill3(X_sag_def(:,j),Y_sag_def(:,j),Z_sag_def(:,j),'r');hold on % first axial plane
% fill3(X_cor_def(:,1),Y_cor_def(:,1),Z_cor_def(:,1),'b');hold on % first sagittal plane
fill3(X_cor_def(:,i),Y_cor_def(:,i),Z_cor_def(:,i),'b');hold on % second sagittal plane
% alpha(.2)

disp('--------- Calculate the source control points ( # N^3 ) -----')
%% Calculate the bounding box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bb = [xmin xmax;ymin ymax;zmin zmax]
% global var_array_v

var_array_v = [var_array1_v;var_array2_v;var_array3_v];

tmp_var1 = var_array_v(:,1);
tmp_var2 = var_array_v(:,2);
tmp_var3 = var_array_v(:,3);

bb = [min(min(tmp_var1(tmp_var1~=-Inf)))-5 max(max(var_array_v(:,1)))+5; ...
      min(min(tmp_var2(tmp_var2~=-Inf)))-5 max(max(var_array_v(:,2)))+5; ...
      min(min(tmp_var3(tmp_var3~=-Inf)))-5 max(max(var_array_v(:,3)))+5];
% X_array = [X_ax_v(:);X_sag_v(:);X_cor_v(:)];
% Y_array = [Y_ax_v(:);Y_sag_v(:);Y_cor_v(:)];
% Z_array = [Z_ax_v(:);Z_sag_v(:);Z_cor_v(:)];
% 
% bb = [min(X_array(:))-25 max(X_array(:))+25;...
%       min(Y_array(:))-25 max(Y_array(:))+25;...
%       min(Z_array(:))-25 max(Z_array(:))+25];
%   
% Create the source control points
nx = 4;
ny = 4;
nz = 3;

l_x = linspace(bb(1,1),bb(1,2),nx);
l_y = linspace(bb(2,1),bb(2,2),ny);
l_z = linspace(bb(3,1),bb(3,2),nz);

global source_control_v
source_control_v = zeros(nx * ny * nz,3);

for i = 1:nx
    for j = 1:ny
        
        tmp =  1:nz;
        s2ind =  tmp + nz*(j-1 + ny*(i-1));
        
        source_control_v(s2ind,1) = repmat(l_x(i),nz,1);
        source_control_v(s2ind,2) = repmat(l_y(j),nz,1);
        source_control_v(s2ind,3) = l_z(tmp);
        
    end
end

% Plot the source control points
for i = 1:nx * ny * nz
    plot3(source_control_v(i,1),source_control_v(i,2),source_control_v(i,3),'k+');hold on
end

disp('--------- Calculate the source control points mesh (tetra_vhedrons) -----')
%% Define the mesh for FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global tetra_v
tetra_v = [];

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
            
            
            tetra_v = [tetra_v;...
                    ind_tmp  ind_tmp2 ind_tmp3 ind_tmp5;...
                    ind_tmp2 ind_tmp3 ind_tmp5 ind_tmp7;...
                    ind_tmp3 ind_tmp7 ind_tmp8 ind_tmp2;...
                    ind_tmp4 ind_tmp7 ind_tmp8 ind_tmp3;...
                    ind_tmp5 ind_tmp7 ind_tmp3 ind_tmp4;...
                    ind_tmp6 ind_tmp7 ind_tmp8 ind_tmp4];
            
            % Just for the plotting
            vari = [source_control_v(ind_tmp,:);source_control_v(ind_tmp2,:);source_control_v(ind_tmp3,:);source_control_v(ind_tmp4,:);...
                    source_control_v(ind_tmp5,:);source_control_v(ind_tmp6,:);source_control_v(ind_tmp7,:);source_control_v(ind_tmp8,:)];
            dt = DelaunayTri(vari);
            
            tetramesh(dt);hold on
            
        end
        
    end
    
end

disp('--------- Convert the computed mesh into DelaunayTri class -----')
%% Convert the computed mesh into DelaunayTri class is gonna help us for computing the
%% vertices or tetra_vhedrons that contain a query of points

global source_tri_v
global list_edges_v

trep = TriRep(tetra_v, source_control_v);
source_tri_v = trep;

alpha(.1)
axis equal

list_edges_v = edges_connected(source_tri_v);

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------- Start Optimization  -----')

disp('--------- Preparing first intersection: axial/sagittal  -----')
size1_int1 = size(var_cell1_v,3);
size2_int1 = size(var_cell1_v,2);
size3_int1 = size(var_cell1_v,1);

global sub_1_int1_v
global sub_2_int1_v
global sub_3_int1_v
global current_tr_int1_v
global c2b_coord_int1_v


[sub_1_int1_v, sub_2_int1_v, sub_3_int1_v] = ind2sub([size1_int1 size2_int1 size3_int1],1:size(var_array1_v,1));
current_tr_int1_v = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int1_v  = cartToBary(source_tri_v,current_tr_int1_v,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri

disp('--------- Preparing second intersection: axial/coronal  -----')
size1_int2 = size(var_cell2_v,3);
size2_int2 = size(var_cell2_v,2);
size3_int2 = size(var_cell2_v,1);

global sub_1_int2_v
global sub_2_int2_v
global sub_3_int2_v
global current_tr_int2_v
global c2b_coord_int2_v


[sub_1_int2_v, sub_2_int2_v, sub_3_int2_v] = ind2sub([size1_int2 size2_int2 size3_int2],1:size(var_array2_v,1));
current_tr_int2_v = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int2_v  = cartToBary(source_tri_v,current_tr_int2_v,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri


disp('--------- Preparing third intersection: sagittal/coronal  -----')
size1_int3 = size(var_cell3_v,3);
size2_int3 = size(var_cell3_v,2);
size3_int3 = size(var_cell3_v,1);

global sub_1_int3_v
global sub_2_int3_v
global sub_3_int3_v
global current_tr_int3_v
global c2b_coord_int3_v


[sub_1_int3_v, sub_2_int3_v, sub_3_int3_v] = ind2sub([size1_int3 size2_int3 size3_int3],1:size(var_array3_v,1));
current_tr_int3_v = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int3_v  = cartToBary(source_tri_v,current_tr_int3_v,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri


lbv = [-20 -20 -20];
ubv = [ 20  20  20];

% lb = repmat(lbv,2*size(source_tri_v.X,1),1) + [source_tri_v.X;source_tri_v.X];
% ub = repmat(ubv,2*size(source_tri_v.X,1),1) + [source_tri_v.X;source_tri_v.X];
lb = [];
ub = [];

tic

preparing_eval;

tim = toc

tic
%% initialization
mesh0 = [source_tri_v.X(:,1:2);source_tri_v.X(:,2:3);source_tri_v.X(:,1) source_tri_v.X(:,3)]; % [source_tri_v.X(:,1:2);source_tri_v.X(:,2:3)]
% 
% mesh0 = [source_tri_v.X;source_tri_v.X;source_tri_v.X];
options = optimset('Display','iter','MaxIter', 40, 'TolFun', .1, 'PlotFcns',@optimplotfval);
[xfinal_tmp fval exitflag output] = fminunc(@myfun_unc_ortho_eval, mesh0, options);
% xfinal = xfinal_tmp;
% 
xfinal(1:size(source_tri_v.X,1),:) = [xfinal_tmp(1:size(source_tri_v.X,1),:) source_tri_v.X(:,3)];
xfinal(1+size(source_tri_v.X,1):2*size(source_tri_v.X,1),:)   = [source_tri_v.X(:,1) xfinal_tmp(1+size(source_tri_v.X,1):2*size(source_tri_v.X,1),:)];
xfinal(1+2*size(source_tri_v.X,1):3*size(source_tri_v.X,1),:) = [ xfinal_tmp(1+2*size(source_tri_v.X,1):3*size(source_tri_v.X,1),1) source_tri_v.X(:,2) xfinal_tmp(1+2*size(source_tri_v.X,1):3*size(source_tri_v.X,1),2)];

%% Define the initial mesh as a vector
% mesh_vector_tmp = [source_tri_v.X;source_tri_v.X]';
% mesh_vector = mesh_vector_tmp(:);
% 
% [xfinal_v, fval] = optimizationWithDE(1, 6*size(source_tri_v.X,1),[],[],[], [], [], [], [], []);
% xfinal = reshape(xfinal_v',2*size(source_tri_v.X,1),3) +  [source_tri_v.X;source_tri_v.X];

% mesh_vector_tmp = [source_tri_v.X;source_tri_v.X;source_tri_v.X]';
% mesh_vector = mesh_vector_tmp(:);

% [xfinal_v, fval] = optimizationWithDE(1, 6 * size(source_tri_v.X,1),[],[],[], [], [], [], [], []);
% [xfinal_v, fval] = optimizationWithDE(1, 9 * size(source_tri_v.X,1),[],[],[], [], [], [], [], []);

% xfinal_t = reshape(xfinal_v', 3 * size(source_tri_v.X,1),3);% +  [source_tri_v.X;source_tri_v.X];

% xfinal(1:size(source_tri_v.X,1),:) = [xfinal_t(1:size(source_tri_v.X,1),:) zeros(size(source_tri_v.X,1),1)] + source_tri_v.X;
% xfinal(size(source_tri_v.X,1)+1:2*size(source_tri_v.X,1),:) = [zeros(size(source_tri_v.X,1),1) xfinal_t(1:size(source_tri_v.X,1),:)] + source_tri_v.X;
% xfinal(2*size(source_tri_v.X,1)+1:3*size(source_tri_v.X,1),:) = [xfinal_t(1:size(source_tri_v.X,1),1) zeros(size(source_tri_v.X,1),1) xfinal_t(1:size(source_tri_v.X,1),2)] + source_tri_v.X;

% target_tri = TriRep(source_tri_v.Triangulation,xfinal);

target_tri_ax  = TriRep(source_tri_v.Triangulation, xfinal(1:size(source_tri_v.X,1),:));
target_tri_sag = TriRep(source_tri_v.Triangulation, xfinal(size(source_tri_v.X,1)+1:2*size(source_tri_v.X,1),:));
target_tri_cor = TriRep(source_tri_v.Triangulation, xfinal(2*size(source_tri_v.X,1)+1:end,:));

% plot_triangulations(source_tri_v,target_tri);
save(save_name,'source_tri_v','target_tri_ax','target_tri_sag','target_tri_cor','xfinal','vol_ax_eval','vol_sag_eval','vol_cor_eval','output');

tim = toc