clc;
close all;

X_ax = [];
Y_ax = [];
Z_ax = [];

X_sag = [];
Y_sag = [];
Z_sag = [];

%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global vol_ax
global vol_sag
global new_im_ax
global new_im_sag

vol_ax = views.axial;
vol_sag = views.sagittal;

new_im_ax = zeros(size(views.axial));
new_im_sag = zeros(size(views.sagittal));

rows = size(views.axial,1);
cols = size(views.axial,2);

total_ax = size(views.axial,3);
total_sag = size(views.sagittal,3);

%figure(1);hold on

global plane_ax
global plane_sag

plane_ax = zeros(total_ax,4);
plane_sag = zeros(total_sag,4);

M = zeros(4,4);
global axial_m
global axial_m1

axial_m = cell(1,total_ax);
axial_m1 = cell(1,total_ax);

for ax=1:size(views.axial,3)
    
    
    M(1:3,4) = views.axial_info{ax}.ImagePositionPatient;
    M(4,4) = 1;
    M(1:3,1) = views.axial_info{ax}.ImageOrientationPatient(1:3).*views.axial_info{ax}.PixelSpacing(1);
    M(1:3,2) = views.axial_info{ax}.ImageOrientationPatient(4:6).*views.axial_info{ax}.PixelSpacing(2);
    
    axial_m{ax} = M;
    axial_m1{ax} = pinv(M);

    [x,y,z] = calculate4corners( M );
    
%     points = [x' y' z'];
%     fill3(points(:,1),points(:,2),points(:,3),'r');hold on
    
    X_ax = [X_ax x'];
    Y_ax = [Y_ax y'];
    Z_ax = [Z_ax z'];
    
    N1 = cross([X_ax(1,ax) Y_ax(1,ax) Z_ax(1,ax)]-[X_ax(2,ax) Y_ax(2,ax) Z_ax(2,ax)],[X_ax(1,ax) Y_ax(1,ax) Z_ax(1,ax)]-[X_ax(3,ax) Y_ax(3,ax) Z_ax(3,ax)]); % normal to the axial (ax)		
	%A1 = [X_ax(4,i) Y_ax(4,i) Z_ax(4,i)];		
    
	plane_ax(ax,:) = [N1 -(N1(1)*X_ax(4,ax) + N1(2)*Y_ax(4,ax) + N1(3)*Z_ax(4,ax))]; % (A,B,C,D) of the plane Ax + By + Cz + d =0
    
end


%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global sag_m
global sag_m1

sag_m = cell(1,total_sag);
sag_m1 = cell(1,total_sag);

for sag=1:size(views.sagittal,3)
    
    
    M(1:3,4) = views.sagittal_info{sag}.ImagePositionPatient;
    M(4,4) = 1;
    M(1:3,1) = views.sagittal_info{sag}.ImageOrientationPatient(1:3).*views.sagittal_info{sag}.PixelSpacing(1);
    M(1:3,2) = views.sagittal_info{sag}.ImageOrientationPatient(4:6).*views.sagittal_info{sag}.PixelSpacing(2);
    
    sag_m{sag} = M;
    sag_m1{sag} = pinv(M);
    
    [x,y,z] = calculate4corners( M );
    
%     points = [x' y' z'];
%     fill3(points(:,1),points(:,2),points(:,3),'b');hold on
    
    X_sag = [X_sag x'];
    Y_sag = [Y_sag y'];
    Z_sag = [Z_sag z'];
    
    N2 = cross([X_sag(1,sag) Y_sag(1,sag) Z_sag(1,sag)]-[X_sag(2,sag) Y_sag(2,sag) Z_sag(2,sag)],[X_sag(1,sag) Y_sag(1,sag) Z_sag(1,sag)]-[X_sag(3,sag) Y_sag(3,sag) Z_sag(3,sag)]);	
	%A2 = [X_sag(4,1) Y_sag(4,1) Z_sag(4,1)];
    
    plane_sag(sag,:) = [N2 -(N2(1)*X_sag(4,sag) + N2(2)*Y_sag(4,sag) + N2(3)*Z_sag(4,sag))]; % (A,B,C,D) of the plane Ax + By + Cz + d =0

end



%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%var = [];

ax = 1:size(views.axial,3);
sag = 1:size(views.sagittal,3);

global t
t = 0:.5:1;
%t = [0 .25 .5 .75 1];

global var_cell
global var_array

var_cell = cell(size(views.axial,3),size(views.sagittal,3),length(t));
var_array = zeros(size(views.axial,3)*size(views.sagittal,3)*length(t),3);

options = optimset('Display','off');

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

        % a = line(V1(:,1),V1(:,2),V1(:,3));
        % line([V1(1,1) V1(1,2) V1(1,3)], [V1(2,1),V1(2,2),V1(2,3)], 'g')

        
        vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];

        for k=1:length(t)
            
            % plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');%,'MarkerSize',i);hold on
            var_cell{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
            
            ind_tmp = sub2ind([size(var_cell,3) size(var_cell,2) size(var_cell,1)],k,j,i);
            %ind_tmp1 = k + length(t)*(j-1 + length(sag)*(i-1));% it is faster than sub2ind

            var_array(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
            var_array(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
            var_array(ind_tmp,3) = V1(1,3) + vd(3)*t(k);

            %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp_v1 = axial_m1{i} * [var_cell{i,j,k} 1]'; % 3D point to 2D point in the frame coordinates
            
            tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            new_i = round(tmp_v(1) + 1); % matlab starts always at 1 
            new_j = round(tmp_v(2) + 1); % matlab starts always at 1
            new_k = i;
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            neig = [views.axial(max(floor(tmp_v(1) + 1),1),max(floor(tmp_v(2) + 1),1),i) views.axial(max(floor(tmp_v(1) + 1),1),min(ceil(tmp_v(2) + 1),cols),i);...
                    views.axial(min(ceil(tmp_v(1) + 1),rows), max(floor(tmp_v(2) + 1),1),i) views.axial(min(ceil(tmp_v(1) + 1),rows), min(ceil(tmp_v(2) + 1),cols),i)];
            new_im_ax(new_i,new_j,new_k) = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
                                                
            
            %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp_v1 = sag_m1{j} * [var_cell{i,j,k} 1]'; % 3D point to 2D point in the frame coordinates
            
            tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            new_i = round(tmp_v(1) + 1); % matlab starts always at 1
            new_j = round(tmp_v(2) + 1); % matlab starts always at 1
            new_k = j;
            
            %% Make sure the indexes are correct
            
            neig = [views.sagittal(max(floor(tmp_v(1) + 1),1),max(floor(tmp_v(2) + 1),1),j) views.sagittal(max(floor(tmp_v(1) + 1),1),min(ceil(tmp_v(2) + 1),cols),j);...
                    views.sagittal(min(ceil(tmp_v(1) + 1),rows), max(floor(tmp_v(2) + 1),1),j) views.sagittal(min(ceil(tmp_v(1) + 1),rows), min(ceil(tmp_v(2) + 1),cols),j)];
            new_im_sag(new_i,new_j,new_k) = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
            
        end
        
       
    end
end

[A_c,b_c,Aeq_c,beq_c] = vert2lcon([var_array(:,1) var_array(:,2) var_array(:,3)]);

%show_results(new_im_ax); % show the lines in the frame coordinates
%show_results(new_im_sag); % show the lines in the frame coordinates

tic
%n_points = 1:length(t);%*size(vol_sag,3);

%[sol,fval] = fminunc(@myfun,var_array(n_points,:));
%[sol,fval] = fminunc(@myfun,ones(size(vol_sag,3)*length(t)*2,3));

%% Prepare the constraints
% cons_Aeq(1:length(n_points),:) = repmat(plane_ax(1,1:3),length(n_points),1);
% cons_beq(1:length(n_points),:) = repmat(plane_ax(1,4),length(n_points),1);
% 
% cons_Aeq(1+length(n_points):2*length(n_points),:) = repmat(plane_sag(1,1:3),length(n_points),1);
% cons_beq(1+length(n_points):2*length(n_points),:) = repmat(plane_sag(1,4),length(n_points),1);

%[sol,fval] = fmincon(@myfun,[var_array(n_points,:)' var_array(n_points,:)']',cons_Aeq,cons_beq);

%% Plot the points

% figure(1);
% for i = 1:length(n_points)
%     plot3(sol(i,1),sol(i,2),sol(i,3),'r*');hold on
%     plot3(var_array(n_points(i),1),var_array(n_points(i),2),var_array(n_points(i),3),'g*');hold on
%     line([sol(i,1),var_array(n_points(i),1)],[sol(i,2),var_array(n_points(i),2)],[sol(i,3),var_array(n_points(i),3)]);hold on
% end

%% Plot the planes

fill3(X_ax(:,1),Y_ax(:,1),Z_ax(:,1),'r');hold on % first axial plane
%fill3(X_ax(:,2),Y_ax(:,2),Z_ax(:,2),'r');hold on % first axial plane
fill3(X_sag(:,1),Y_sag(:,1),Z_sag(:,1),'b');hold on % first sagittal plane
%fill3(X_sag(:,2),Y_sag(:,2),Z_sag(:,2),'b');hold on % second sagittal plane    
alpha(.1)
axis equal

time = toc

%% Define the mesh for FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure; % plot the mesh once the intersection points are defined

% for i=1:length(ax)-1
%     for j=1:length(sag)-1
% 
%         for k=1:length(t)-1
%            
%             vari = [var_cell{i,j,k};var_cell{i,j,k+1};var_cell{i+1,j,k};var_cell{i+1,j+1,k};var_cell{i,j+1,k};var_cell{i+1,j+1,k+1};var_cell{i,j+1,k+1};var_cell{i+1,j,k+1}];
%             ind_tmp = sub2ind([size(var_cell,3) size(var_cell,2) size(var_cell,1)],k,j,i);

%             tmp = calculate_tetrahedron(i,j,k,t,sag);
%             tetra = [tetra;tmp];

%             %dt = DelaunayTri(vari);
%             % ind2sub(size(var_cell),i,j,k);
%             %vert_mat(k,j,i) = (i-1)*(length(ax)-1) + (k-1)+ (j-1) + (k-1)*(length(sag)-1);
%             %plot3(var_cell{i,j,k}(1), var_cell{i,j,k}(2), var_cell{i,j,k}(3),'g*');hold on
%             %tetramesh(dt,vari,'facealpha',.3);hold on
%                            
%         end
%                 
%     end
% end

%% Calculate the bounding box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bb = [xmin xmax;ymin ymax;zmin zmax]

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

for i=1:nx
    for j=1:ny
        for k=1:nz
            s2ind =  k + nz*(j-1 + ny*(i-1));
            source_control(s2ind,:) = [l_x(i) l_y(j) l_z(k)];
        end
    end
end


for i = 1:nx * ny * nz
    plot3(source_control(i,1),source_control(i,2),source_control(i,3),'k+');hold on
end






