function f = myfun(list_x)

global var_cell
global sag_m1
global axial_m1
global vol_ax
global vol_sag
global t

rows = size(vol_ax,1);
cols = size(vol_ax,2);
f = 0;

for i = 1:2%size(vol_ax,3)
   for j = 1:size(vol_sag,3)
       for k = 1:length(t)
           
           %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ind_tmp = sub2ind([size(var_cell,3) size(var_cell,2) size(var_cell,1)],k,j,i);
            %ind_tmp = ind_tmp - length(t)*size(vol_sag,3);
            %ind_tmp = length(t)*(j-1) + k + (i-1)*cols*length(t);

            tmp_v1 = axial_m1{i} * [list_x(ind_tmp,:) 1]'; % 3D point to 2D point in the frame coordinates
            
            tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            %% Make sure the indexes are correct and use bilinear interpolation
            
            neig = [vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(floor(tmp_v(2) + 1),1),cols),i) vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(ceil(tmp_v(2) + 1),1),cols),i);...
                    vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(floor(tmp_v(2) + 1),1),cols),i) vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(ceil(tmp_v(2) + 1),1),cols),i)];
            new_im_ax = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
                                                
            
            %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp_v1 = sag_m1{j} * [var_cell{i,j,k} 1]'; % 3D point to 2D point in the frame coordinates 
            
            tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
            
            %% Make sure the indexes are correct
            
            neig = [vol_sag(min(max(floor(tmp_v(1) + 1),1),rows),min(max(floor(tmp_v(2) + 1),1),cols),j) vol_sag(min(max(floor(tmp_v(1) + 1),1),rows),min(max(ceil(tmp_v(2) + 1),1),cols),j);...
                    vol_sag(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(floor(tmp_v(2) + 1),1),cols),j) vol_sag(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(ceil(tmp_v(2) + 1),1),cols),j)];
            new_im_sag = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
            
            %% Function
            f = f + (new_im_ax -  new_im_sag)^2;
           
       end
   end
end
