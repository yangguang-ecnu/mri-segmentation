function [F J] = myfun3(mesh0_X, in, out)

tri1 = TriRep(in.source_tri.Triangulation,mesh0_X(1:in.size_source,:)); % define the new mesh for the axial
tri2 = TriRep(in.source_tri.Triangulation,mesh0_X(in.size_source+1:end,:)); % define the new mesh for the sagittal


b2c_ncoord1 = baryToCart(tri1, in.current_tr, in.c2b_coord);
b2c_ncoord2 = baryToCart(tri2, in.current_tr, in.c2b_coord);

 for i = 1 : in.n  
     
           %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tmp_v1 = in.axial_m1{in.sub_3(i)} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
            
            %tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            %% Make sure the indexes are correct and use bilinear interpolation
            fl = floor(tmp_v1(2) + 1);
            fl2 = floor(tmp_v1(1) + 1);
            cl = ceil(tmp_v1(1) + 1);
            cl2 = ceil(tmp_v1(2) + 1);
            
            neig = [in.vol_ax(min(max(fl,1),in.rows),min(max(fl2,1),in.cols),in.sub_3(i))   in.vol_ax(min(max(fl,1),in.rows),min(max(cl,1),in.cols),in.sub_3(i));...
                    in.vol_ax(min(max(cl2,1),in.rows), min(max(fl2,1),in.cols),in.sub_3(i)) in.vol_ax(min(max(cl2,1),in.rows), min(max(cl,1),in.cols),in.sub_3(i))];
            new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1), double(neig));
                                                
            grad_s1 = [sum(neig(:,1)) - sum(neig(:,2)) sum(neig(1,:)) - sum(neig(2,:))];%[Gx Gy]
            
            %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tmp_v1 = in.sag_m1{in.sub_2(i)} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
            
            %tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            %% Make sure the indexes are correct
            fl = floor(tmp_v1(2) + 1);
            fl2 = floor(tmp_v1(1) + 1);
            cl = ceil(tmp_v1(1) + 1);
            cl2 = ceil(tmp_v1(2) + 1);
            
            neig = [in.vol_sag(min(max(fl,1),in.rows),min(max(fl2,1),in.cols),in.sub_2(i))   in.vol_sag(min(max(fl,1),in.rows),min(max(cl,1),in.cols),in.sub_2(i));...
                    in.vol_sag(min(max(cl2,1),in.rows), min(max(fl2,1),in.cols),in.sub_2(i)) in.vol_sag(min(max(cl2,1),in.rows), min(max(cl,1),in.cols),in.sub_2(i))];
                    
            new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1), double(neig));
            
            grad_s2 = [sum(neig(:,1)) - sum(neig(:,2)) sum(neig(1,:)) - sum(neig(2,:))];%[Gx Gy]
            
            %% Jacobian of warp function qi = w(pi,X), pi is a given data point, and X is the mesh
            in.w1(1,(in.source_tri.Triangulation(in.current_tr(i),:)-1)*3 + 1) = in.c2b_coord(i,:);
            in.w1(2,:) = circshift(in.w1(1,:)',1);
            in.w1(3,:) = circshift(in.w1(1,:)',2);
            
            %% Function

           out.F(i) = new_im_ax - new_im_sag;
           in.J1(i,:) = grad_s1 * in.axial_m1{in.sub_3(i)}(1:2,1:3) * in.w1;  
           in.J2(i,:) = grad_s2 * in.sag_m1{in.sub_2(i)}(1:2,1:3)   * in.w1;
 end
 
 out.J(1:in.n,:) = [in.J1 in.J2];


 
 % Set the matrix to zero
 in.w1 = in.w1 - in.w1;
 
 
 
 for i = 1:in.size_source
     
     in.F(i + in.n) = in.lambda*(norm( in.source_tri.X(i,:) - mean(in.source_tri.X(in.list_edges{i},:)) ));
     mean_tmp = mean(in.source_tri.X(in.list_edges{i},:));
     
     in.J(i + in.n,(i-1)*3 + 1:(i-1)*3 + 3) = [2*in.lambda*(in.source_tri.X(i,1) - mean_tmp(1) ) 2*in.lambda*(in.source_tri.X(i,2) - mean_tmp(2) ) 2*in.lambda*(in.source_tri.X(i,3) - mean_tmp(3) )];
     in.J(i + in.n,(in.list_edges{i}-1)*3 + 1:(in.list_edges{i}-1)*3 + 3) = [in.scalar(i) in.scalar(i) in.scalar(i)];
 end
 
 J = out.J;
 F = out.F;