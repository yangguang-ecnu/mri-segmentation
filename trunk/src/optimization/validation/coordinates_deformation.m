function [out] = coordinates_deformation(t, st, n_x, n_y, n_z, sigma)

%% Calculate the intersection between the different planes

options = optimset('Display','off');

ax  = st.slices_ax;
sag = st.slices_sag;
cor = st.slices_cor;


var_cell1_v  = cell(length(st.slices_ax),length(st.slices_sag),length(t));
var_array1_v = zeros(length(st.slices_ax)*length(st.slices_sag)*length(t),3);

disp('--------- Calculate the s1 x s2 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(ax)
    
    for j=1:length(sag)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([st.X_ax(:,i)  st.Y_ax(:,i)  st.Z_ax(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([st.X_sag(:,j) st.Y_sag(:,j) st.Z_sag(:,j)]); %% constraints of the sagittal plane
        
        normal = cross([st.X_sag(1,j) st.Y_sag(1,j) st.Z_sag(1,j)]-[st.X_sag(2,j) st.Y_sag(2,j) st.Z_sag(2,j)],[st.X_sag(1,j) st.Y_sag(1,j) st.Z_sag(1,j)]-[st.X_sag(3,j) st.Y_sag(3,j) st.Z_sag(3,j)]); % direction of the intersection line
        
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
                
%                 plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'r+');hold on;%,'MarkerSize',i);hold on
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


var_cell2_v = cell(length(st.slices_ax),length(st.slices_cor),length(t));
var_array2_v = zeros(length(st.slices_ax)*length(st.slices_cor)*length(t),3);

options = optimset('Display','off');


for i=1:length(ax)
    
    for j=1:length(cor)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([st.X_ax(:,i) st.Y_ax(:,i) st.Z_ax(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([st.X_cor(:,j) st.Y_cor(:,j) st.Z_cor(:,j)]); %% constraints of the sagittal plane
        
        if i==1 && j==1
            normal = cross([st.X_cor(1,j) st.Y_cor(1,j) st.Z_cor(1,j)]-[st.X_cor(2,j) st.Y_cor(2,j) st.Z_cor(2,j)],[st.X_cor(1,j) st.Y_cor(1,j) st.Z_cor(1,j)]-[st.X_cor(3,j) st.Y_cor(3,j) st.Z_cor(3,j)]); % direction of the intersection line
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
                
%                 plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');hold on%,'MarkerSize',i);hold on
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


% fill3(st.X_ax(:,1),st.Y_ax(:,1),st.Z_ax(:,1),'r');hold on % first axial plane
% fill3(st.X_ax(:,total_ax),st.Y_ax(:,total_ax),st.Z_ax(:,total_ax),'r');hold on % first axial plane
% fill3(st.X_cor(:,1),st.Y_cor(:,1),st.Z_cor(:,1),'b');hold on % first sagittal plane
% fill3(st.X_cor(:,total_cor),st.Y_cor(:,total_cor),st.Z_cor(:,total_cor),'b');hold on % second sagittal plane
% alpha(.2)
disp('--------- Calculate the s2 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


var_cell3_v  = cell(length(st.slices_cor),length(st.slices_sag),length(t));
var_array3_v = zeros(length(st.slices_sag)*length(st.slices_cor)*length(t),3);

options = optimset('Display','off');

for i=1:length(cor)
    
    for j=1:length(sag)
        
        
        
        [A1,b1,Aeq1,beq1] = vert2lcon([st.X_cor(:,i) st.Y_cor(:,i) st.Z_cor(:,i)]); %% constraints of the axial plane
        [A2,b2,Aeq2,beq2] = vert2lcon([st.X_sag(:,j) st.Y_sag(:,j) st.Z_sag(:,j)]); %% constraints of the sagittal plane
        
        if i==1 && j==1
            normal = cross([st.X_sag(1,j) st.Y_sag(1,j) st.Z_sag(1,j)]-[st.X_sag(2,j) st.Y_sag(2,j) st.Z_sag(2,j)],[st.X_sag(1,j) st.Y_sag(1,j) st.Z_sag(1,j)]-[st.X_sag(3,j) st.Y_sag(3,j) st.Z_sag(3,j)]); % direction of the intersection line
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
                
%                 plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'b+'); hold on;
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


% c' {ijk} = c{ijk} + sigma * randn , sigma = deformable amplitude

X_array = [st.X_ax_t(:);st.X_sag_t(:);st.X_cor_t(:)];
Y_array = [st.Y_ax_t(:);st.Y_sag_t(:);st.Y_cor_t(:)];
Z_array = [st.Z_ax_t(:);st.Z_sag_t(:);st.Z_cor_t(:)];

bb = [min(X_array(:))-5 max(X_array(:))+5; ...
      min(Y_array(:))-5 max(Y_array(:))+5; ...
      min(Z_array(:))-5 max(Z_array(:))+5];

% Create the source control points

l_x = linspace(bb(1,1),bb(1,2),n_x);
l_y = linspace(bb(2,1),bb(2,2),n_y);
l_z = linspace(bb(3,1),bb(3,2),n_z);


control_points = zeros(n_x * n_y * n_z, 3);

for i = 1:n_x
    for j = 1:n_y
        
        tmp =  1:n_z;
        s2ind =  tmp + n_z*(j-1 + n_y*(i-1));
        
        control_points(s2ind,1) = repmat(l_x(i),n_z,1);
        control_points(s2ind,2) = repmat(l_y(j),n_z,1);
        control_points(s2ind,3) = l_z(tmp);
        
    end
end

% Plot the source control points
% for i = 1:n_x * n_y * n_z
%     plot3(control_points(i,1),control_points(i,2),control_points(i,3),'k+');hold on
% end

disp('--------- Calculate the source control points mesh (tetra_vhedrons) -----')
%% Define the mesh for FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tetra_c = [];

% figure;
for i=1:n_x-1
    
    for j=1:n_y-1
        
        for k=1:n_z-1
            
            ind_tmp =  k + n_z*(j-1 + n_y*(i-1));
            
            ind_tmp2 =  k+1 + n_z*(j-1 + n_y*(i-1));
            ind_tmp3 =  k + n_z*(j-1 + n_y*(i));
            ind_tmp4 =  k + n_z*(j + n_y*(i));
            ind_tmp5 =  k + n_z*(j + n_y*(i-1));
            ind_tmp6 =  k+1 + n_z*(j + n_y*(i));
            ind_tmp7 =  k+1 + n_z*(j + n_y*(i-1));
            ind_tmp8 =  k+1 + n_z*(j-1 + n_y*(i));
            
            
            tetra_c = [tetra_c;...
                ind_tmp  ind_tmp2 ind_tmp3 ind_tmp5;...
                ind_tmp2 ind_tmp3 ind_tmp5 ind_tmp7;...
                ind_tmp3 ind_tmp7 ind_tmp8 ind_tmp2;...
                ind_tmp4 ind_tmp7 ind_tmp8 ind_tmp3;...
                ind_tmp5 ind_tmp7 ind_tmp3 ind_tmp4;...
                ind_tmp6 ind_tmp7 ind_tmp8 ind_tmp4];
            
            % Just for the plotting
%             vari = [control_points(ind_tmp,:);control_points(ind_tmp2,:);control_points(ind_tmp3,:);control_points(ind_tmp4,:);...
%                 control_points(ind_tmp5,:);control_points(ind_tmp6,:);control_points(ind_tmp7,:);control_points(ind_tmp8,:)];
%             dt = DelaunayTri(vari);
%             
%             tetramesh(dt);hold on
            
        end
        
    end
    
end
% alpha(.1)
% axis equal

disp('--------- Convert the computed mesh into DelaunayTri class -----')
%% Convert the computed mesh into DelaunayTri class is gonna help us for computing the
%% vertices or tetra_vhedrons that contain a query of points

trep = TriRep(tetra_c, control_points);
control_tri = trep;

% Deformed grid
 % it is the deformation amplitude
deform_points_ax  = control_points(:,1:2) + sigma .* randn(size(control_points(:,1:2)));
deform_points_sag = control_points(:,2:3) + sigma .* randn(size(control_points(:,1:2)));
deform_points_cor = [control_points(:,1) control_points(:,3)] + sigma .* randn(size(control_points(:,1:2)));

deform_points_ax3  = [deform_points_ax control_tri.X(:,3)];
deform_points_sag3 = [control_tri.X(:,1) deform_points_sag];
deform_points_cor3 = [deform_points_cor(:,1) control_tri.X(:,2) deform_points_cor(:,2)];

deform_tri_ax  = TriRep(tetra_c, deform_points_ax3);
deform_tri_sag = TriRep(tetra_c, deform_points_sag3);
deform_tri_cor = TriRep(tetra_c, deform_points_cor3);


%% Plot the two grids
% figure;
% subplot(121);tetramesh(control_tri);hold on
% subplot(122);tetramesh(deform_tri);hold on
% 
% alpha(.1)
% axis equal

disp('--------- Compute the new images -----')
%% Compute the new images


current_tr1 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
current_tr2 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]);
current_tr3 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]);

% New data positions with respect the new computed meshes
b1 = cartToBary(control_tri,  current_tr1,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc1_ax  = baryToCart(deform_tri_ax,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
xc1_sag = baryToCart(deform_tri_sag, current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )

b2 = cartToBary(control_tri,  current_tr2,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc2_ax  = baryToCart(deform_tri_ax, current_tr2, b2); % cartesian coordinates of the points wrt target tri ( p' )
xc2_cor = baryToCart(deform_tri_cor, current_tr2, b2); % cartesian coordinates of the points wrt target tri ( p' )

b3 = cartToBary(control_tri,  current_tr3,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc3_sag = baryToCart(deform_tri_sag, current_tr3, b3); % cartesian coordinates of the points wrt target tri ( p' )
xc3_cor = baryToCart(deform_tri_cor, current_tr3, b3); % cartesian coordinates of the points wrt target tri ( p' )

out.xc1_ax  = xc1_ax;
out.xc1_sag = xc1_sag;
out.xc2_ax  = xc2_ax;
out.xc2_cor = xc2_cor;
out.xc3_sag = xc3_sag;
out.xc3_cor = xc3_cor;

out.deform_tri_ax  = deform_tri_ax;
out.deform_tri_sag = deform_tri_sag;
out.deform_tri_cor = deform_tri_cor;

out.var_array3 = var_array3_v;
out.var_array2 = var_array2_v;
out.var_array1 = var_array1_v;

out.var_cell3 = var_cell3_v;
out.var_cell2 = var_cell2_v;
out.var_cell1 = var_cell1_v;

out.control_tri = control_tri;