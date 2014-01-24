function F = myfun2(mesh0_X)

global source_tri
%global var_cell1
global var_array1
global sag_m1
global axial_m1
global vol_ax
global vol_sag
%global sub_1
global sub_2
global sub_3
%global t
global current_tr
global c2b_coord
global list_edges

rows = size(vol_ax,1);
cols = size(vol_ax,2);


%tri = TriRep(source_tri.Triangulation,mesh0_X);

tri1 = TriRep(source_tri.Triangulation,mesh0_X(1:size(source_tri.X,1),:)); % define the new mesh for the axial
tri2 = TriRep(source_tri.Triangulation,mesh0_X(size(source_tri.X,1)+1:end,:)); % define the new mesh for the sagittal
% 
% tri1 = TriRep(source_tri.Triangulation,source_tri.X + mesh0_X(1:size(source_tri.X,1),:)); % define the new mesh for the axial
% tri2 = TriRep(source_tri.Triangulation,source_tri.X + mesh0_X(size(source_tri.X,1)+1:end,:)); % define the new mesh for the sagittal


F = zeros(size(var_array1,1),1);
%J = zeros(size(var_array1,1),size(source_tri.X,1));

% w1 = zeros(3,3*size(source_tri.X,1));
% w2 = zeros(3,3*size(source_tri.X,1));

b2c_ncoord1 = baryToCart(tri1, current_tr, c2b_coord);
b2c_ncoord2 = baryToCart(tri2, current_tr, c2b_coord);

n = size(var_array1,1);

 for i = 1:n  
     
           %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tmp_v1 = axial_m1{sub_3(i)} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
            
            tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            neig = [vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(floor(tmp_v(2) + 1),1),cols),sub_3(i)) vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(ceil(tmp_v(2) + 1),1),cols),sub_3(i));...
                    vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(floor(tmp_v(2) + 1),1),cols),sub_3(i)) vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(ceil(tmp_v(2) + 1),1),cols),sub_3(i))];
            new_im_ax = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
                                                
            %grad_s1 = [sum(neig(:,1)) - sum(neig(:,2)) sum(neig(1,:)) - sum(neig(2,:))];%[Gx Gy]
            
            %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tmp_v1 = sag_m1{sub_2(i)} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
            
            %tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            %% Make sure the indexes are correct
            
            neig = [vol_sag(min(max(floor(tmp_v(1) + 1),1),rows),min(max(floor(tmp_v(2) + 1),1),cols),sub_2(i)) vol_sag(min(max(floor(tmp_v(1) + 1),1),rows),min(max(ceil(tmp_v(2) + 1),1),cols),sub_2(i));...
                    vol_sag(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(floor(tmp_v(2) + 1),1),cols),sub_2(i)) vol_sag(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(ceil(tmp_v(2) + 1),1),cols),sub_2(i))];
            new_im_sag = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
            
%             grad_s2 = [sum(neig(:,1)) - sum(neig(:,2)) sum(neig(1,:)) - sum(neig(2,:))];%[Gx Gy]
%             
%             %% Jacobian of warp function qi = w(pi,X), pi is a given data point, and X is the mesh
%             w1(1,(source_tri.Triangulation(current_tr(i),:)-1)*3 + 1) = c2b_coord(i,:);
%             w1(2,:) = circshift(w1(1,:)',1);
%             w1(3,:) = circshift(w1(1,:)',2);
            
%            w2 = w1;
            
%             w2(1,source_tri.Triangulation(current_tr(i),:) + [0 3 6 9]) = c2b_coord(i,:);
%             w2(2,:) = circshift(w1(1,:)',1);
%             w2(3,:) = circshift(w1(1,:)',2);
            
            %% Function
%            size(grad_s1)
%            size(axial_m1{sub_3(i)}(1:2,1:3))
%            size(w1)
             F(i) = new_im_ax - new_im_sag;
%            J1(i,:) = grad_s1 * axial_m1{sub_3(i)}(1:2,1:3) * w1;  
%            J2(i,:) = grad_s2 * sag_m1{sub_2(i)}(1:2,1:3)   * w1;
 end
 
% output.J = [J1 J2];

%  lambda = .5;
%  for i = 1:size(source_tri.X,1)
%      
%      F(i + n) = lambda*(norm( source_tri.X(i,:) - mean(source_tri.X(list_edges{i},:)) ));
%      
%  end

 
 
 