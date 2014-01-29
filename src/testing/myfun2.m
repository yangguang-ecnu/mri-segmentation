function F = myfun2(mesh0_X,in,out)

tri1 = TriRep(in.source_tri.Triangulation,mesh0_X(1:size(in.source_tri.X,1),:)); % define the new mesh for the axial
tri2 = TriRep(in.source_tri.Triangulation,mesh0_X(size(in.source_tri.X,1)+1:end,:)); % define the new mesh for the sagittal


b2c_ncoord1 = baryToCart(tri1, in.current_tr, in.c2b_coord);
b2c_ncoord2 = baryToCart(tri2, in.current_tr, in.c2b_coord);

incr = in.n + in.size_source;

 for i = 1:in.n  
     
           %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tmp_v1 = in.axial_m1{in.sub_3(i)} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
            
            tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            neig = [in.vol_ax(min(max(floor(tmp_v(1) + 1),1),in.rows),min(max(floor(tmp_v(2) + 1),1),in.cols),in.sub_3(i)) in.vol_ax(min(max(floor(tmp_v(1) + 1),1),in.rows),min(max(ceil(tmp_v(2) + 1),1),in.cols),in.sub_3(i));...
                    in.vol_ax(min(max(ceil(tmp_v(1) + 1),1),in.rows), min(max(floor(tmp_v(2) + 1),1),in.cols),in.sub_3(i)) in.vol_ax(min(max(ceil(tmp_v(1) + 1),1),in.rows), min(max(ceil(tmp_v(2) + 1),1),in.cols),in.sub_3(i))];
            new_im_ax = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
                                                

            
            %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tmp_v1 = in.sag_m1{in.sub_2(i)} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
            
            tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)];
            %% Make sure the indexes are correct
            
            neig = [in.vol_sag(min(max(floor(tmp_v(2) + 1),1),in.rows),min(max(floor(tmp_v(1) + 1),1),in.cols),in.sub_2(i)) in.vol_sag(min(max(floor(tmp_v(2) + 1),1),in.rows),min(max(ceil(tmp_v(1) + 1),1),in.cols),in.sub_2(i));...
                    in.vol_sag(min(max(ceil(tmp_v(2) + 1),1),in.rows), min(max(floor(tmp_v(1) + 1),1),in.cols),in.sub_2(i)) in.vol_sag(min(max(ceil(tmp_v(2) + 1),1),in.rows), min(max(ceil(tmp_v(1) + 1),1),in.cols),in.sub_2(i))];
            new_im_sag = bilinear_interpolation(tmp_v(2),tmp_v(1),double(neig));
            

            %% Function
            out.F(i) = new_im_ax - new_im_sag;

 end
 

 for i = 1:size(in.source_tri.X,1)
     
     out.F(i + in.n) = in.lambda*(sum( (tri1.X(i,:) - mean(tri2.X(in.list_edges{i},:))).^ 2 ));
     out.F(i + incr) = in.lambda*(sum( (tri2.X(i,:) - mean(tri2.X(in.list_edges{i},:))).^ 2 ));
 end

 F = out.F;
 
 