function [out] = coordinates_deformation_rot(t,st,control_tri,def)

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


%% Compute the new images
current_tr1 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
current_tr2 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]);
current_tr3 = tsearchn(control_tri.X,control_tri.Triangulation,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]);

% New data positions with respect the new computed meshes
b1 = cartToBary(control_tri,  current_tr1,[var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc1_ax  = baryToCart(def.deform_tri_ax,  current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )
xc1_sag = baryToCart(def.deform_tri_sag, current_tr1, b1); % cartesian coordinates of the points wrt target tri ( p' )

b2 = cartToBary(control_tri,  current_tr2,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc2_ax  = baryToCart(def.deform_tri_ax, current_tr2, b2); % cartesian coordinates of the points wrt target tri ( p' )
xc2_cor = baryToCart(def.deform_tri_cor, current_tr2, b2); % cartesian coordinates of the points wrt target tri ( p' )

b3 = cartToBary(control_tri,  current_tr3,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % barycentric coordinates of the points wrt source tri
xc3_sag = baryToCart(def.deform_tri_sag, current_tr3, b3); % cartesian coordinates of the points wrt target tri ( p' )
xc3_cor = baryToCart(def.deform_tri_cor, current_tr3, b3); % cartesian coordinates of the points wrt target tri ( p' )

out.xc1_ax  = xc1_ax;
out.xc1_sag = xc1_sag;
out.xc2_ax  = xc2_ax;
out.xc2_cor = xc2_cor;
out.xc3_sag = xc3_sag;
out.xc3_cor = xc3_cor;

out.var_array3 = var_array3_v;
out.var_array2 = var_array2_v;
out.var_array1 = var_array1_v;

out.var_cell3 = var_cell3_v;
out.var_cell2 = var_cell2_v;
out.var_cell1 = var_cell1_v;

