function [var_cell, var_array] = calculate_intersections(X1, Y1, Z1, X2, Y2, Z2, t, size1, size2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_cell  = cell(size1,size2,length(t));
var_array = zeros(size1*size2*length(t),3);

options = optimset('Display','off');

for i=1:size(X1,2)
    
    for j=1:size(X2,2)
        
        [A1, b1, Aeq1, beq1] = vert2lcon([X1(:,i) Y1(:,i) Z1(:,i)]); %% constraints of the axial plane
        [A2, b2, Aeq2, beq2] = vert2lcon([X2(:,j) Y2(:,j) Z2(:,j)]); %% constraints of the sagittal plane
        
        if i==1 && j==1
            normal = cross([X2(1,j) Y2(1,j) Z2(1,j)]-[X2(2,j) Y2(2,j) Z2(2,j)],[X2(1,j) Y2(1,j) Z2(1,j)]-[X2(3,j) Y2(3,j) Z2(3,j)]); % direction of the intersection line
        end
        
        %% Concatenate both constraints
        A = [A1;A2];
        b = [b1;b2];
        Aeq = [Aeq1;Aeq2];
        beq = [beq1;beq2];
        
        [x0,~,~,~,~] = linprog(normal,A,b,Aeq,beq,[],[],[],options); % use linear programming to determine one solution
        
        [V1,~,~] = qlcon2vert(x0, A, b, Aeq, beq); % given the constraints and  a solution, determine the vertices, the intersection points !!
        
%         plot3(V1(:,1),V1(:,2),V1(:,3),'g*');hold on
         
        vd = - [V1(1,1),V1(1,2),V1(1,3)] + [V1(2,1),V1(2,2),V1(2,3)];
        
        for k=1:length(t)
            
%             plot3(V1(1,1) + vd(1)*t(k),V1(1,2) + vd(2)*t(k),V1(1,3) + vd(3)*t(k),'g+');
            var_cell{i,j,k} = [V1(1,1) + vd(1)*t(k) V1(1,2) + vd(2)*t(k) V1(1,3) + vd(3)*t(k)]; %% the intersection points
            
            ind_tmp = sub2ind([size(var_cell,3) size(var_cell,2) size(var_cell,1)],k,j,i);
            
            var_array(ind_tmp,1) = V1(1,1) + vd(1)*t(k);
            var_array(ind_tmp,2) = V1(1,2) + vd(2)*t(k);
            var_array(ind_tmp,3) = V1(1,3) + vd(3)*t(k);
            
        end
        
        
    end
end