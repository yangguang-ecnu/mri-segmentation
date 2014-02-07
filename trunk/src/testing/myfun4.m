function [F J] = myfun4(mesh0_X, in)

tri1 = TriRep(in.source_tri.Triangulation,mesh0_X(1:in.size_source,:)); % define the new mesh for the axial
tri2 = TriRep(in.source_tri.Triangulation,mesh0_X(in.size_source+1:end,:)); % define the new mesh for the sagittal


b2c_ncoord1 = baryToCart(tri1, in.current_tr, in.c2b_coord);
b2c_ncoord2 = baryToCart(tri2, in.current_tr, in.c2b_coord);

incr_row = in.n + in.size_source;
incr_col = 3 * in.size_source;

sub_3 = in.sub_3;
sub_2 = in.sub_2;
rows = in.rows;
cols = in.cols;
vol_ax = in.vol_ax;
vol_sag = in.vol_sag;
F = zeros(in.n + 2 * in.size_source, 1);

% J1 = zeros(in.n, 3 * in.size_source);
% J2 = zeros(in.n, 3 * in.size_source);
source_tri = in.source_tri;
current_tr = in.current_tr;
c2b_coord = in.c2b_coord;
axi_m1_par = in.axial_m1_par;
sag_m1_par = in.sag_m1_par;
%size_w1 = size(in.w1);
neigax = in.neig_ax;
neigsg = in.neig_sag;

for i = 1 : in.n
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_ax = axi_m1_par{i} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
    %% Make sure the indexes are correct and use bilinear interpolation
    
    fl1ax = floor(tmp_v1_ax(1) + 1);
    fl2ax = floor(tmp_v1_ax(2) + 1);
    cl1ax = ceil(tmp_v1_ax(1) + 1);
    cl2ax = ceil(tmp_v1_ax(2) + 1);
    
    min_max_r1ax = min(max(fl2ax,1),rows);
    min_max_r2ax = min(max(cl2ax,1),rows);
    min_max_c1ax = min(max(fl1ax,1),cols);
    min_max_c2ax = min(max(cl1ax,1),cols);
    
    neigax = [vol_ax(min_max_r1ax, min_max_c1ax, sub_3(i)) vol_ax(min_max_r1ax, min_max_c2ax, sub_3(i));...
              vol_ax(min_max_r2ax, min_max_c1ax, sub_3(i)) vol_ax(min_max_r2ax, min_max_c2ax, sub_3(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax));
    
%    grad_s1 = [sum(neigax(:,2)) - sum(neigax(:,1)) sum(neigax(1,:)) - sum(neigax(2,:))];
    
    %neigax = block([min_max_r1ax min_max_c1ax],[5 5],vol_ax(:,:,sub_3(i)));
%     
%     new_im_ax = neigax{sub_3(i)}(min_max_r1ax, min_max_c1ax);
%     
%     [grad_s1_x grad_s1_y] = gradient(neigax);%[sum(neigax(:,2)) - sum(neigax(:,1)) sum(neigax(1,:)) - sum(neigax(2,:))];%[Gx Gy]
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_sg = sag_m1_par{i} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1{sub_2(i)}
    
    %% Make sure the indexes are correct
    
    fl1sg = floor(tmp_v1_sg(1) + 1);
    fl2sg = floor(tmp_v1_sg(2) + 1);
    cl1sg = ceil(tmp_v1_sg(1) + 1);
    cl2sg = ceil(tmp_v1_sg(2) + 1);
    
    min_max_r1sg = min(max(fl2sg,1),rows);
    min_max_r2sg = min(max(cl2sg,1),rows);
    min_max_c1sg = min(max(fl1sg,1),cols);
    min_max_c2sg = min(max(cl1sg,1),cols);
    
    neigsg = [vol_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) vol_sag(min_max_r1sg, min_max_c2sg,sub_2(i));...
              vol_sag(min_max_r2sg, min_max_c1sg, sub_2(i)) vol_sag(min_max_r2sg, min_max_c2sg,sub_2(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg));
    
%    grad_s2 = [sum(neigsg(:,2)) - sum(neigsg(:,1)) sum(neigsg(1,:)) - sum(neigsg(2,:))];
%     neigsg = block([min_max_r1sg min_max_c1sg],[5 5],vol_sag(:,:,sub_2(i)));
    
%     new_im_sag = neigsg{sub_2(i)}(min_max_r1sg, min_max_c1sg);
%     
%     [grad_s2_x grad_s2_y] = gradient(neigsg);
    
    %% Jacobian of warp function qi = w(pi,X), pi is a given data point, and X is the mesh
%     w1 = zeros(size_w1);
%     
%     w1(1,(source_tri.Triangulation(current_tr(i),:)-1)*3 + 1) = c2b_coord(i,:); % source_tri.Triangulation(current_tr(i),:)-1)
%     w1(2,:) = circshift(w1(1,:)',1);
%     w1(3,:) = circshift(w1(1,:)',2);
%     
%     w1(4, w1(2,:) | w1(3,:) | w1(1,:)) = 1;
%     
    w  = zeros(3, 2 * incr_col);
    w2 = zeros(3, 2 * incr_col);
    
    ind_w = (source_tri.Triangulation(current_tr(i),:)-1)*3 + 1;
    
    w(1,ind_w )            = c2b_coord(i,:);
    w2(1,ind_w + incr_col) = c2b_coord(i,:);
    
    w(2,:) = circshift(w(1,:)',1);
    w(3,:) = circshift(w(1,:)',2);
    w2(2,:) = circshift(w2(1,:)',1);
    w2(3,:) = circshift(w2(1,:)',2);
    
%     w( 4, w(2,:)  | w(3,:)  | w(1,:)) = 1;
%     w2(4, w2(2,:) | w2(3,:) | w2(1,:)) = 1;
    %% Function
    
    F(i) = new_im_ax - new_im_sag;
    
%     J1(i,:) = [grad_s1 0] * axi_m1_par{i} * w1; % axial_m1{sub_3(i)}(1:2,:)
%     J2(i,:) = [grad_s2 0] * sag_m1_par{i} * w1; % sag_m1{sub_2(i)}(1:2,:)
    
    grad_ax  = [in.gradx_ax( min_max_r1ax, min_max_c1ax, sub_3(i)) in.grady_ax( min_max_r1ax, min_max_c1ax, sub_3(i))];
    grad_sag = [in.gradx_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) in.grady_sag(min_max_r1sg, min_max_c1sg, sub_2(i))];
    
    Js(i,:) = (grad_ax * axi_m1_par{i}(1:2,1:3) * w) - (grad_sag * sag_m1_par{i}(1:2,1:3) * w2);
end

% Js = [J1 J2];

for i = 1:in.size_source
    
    mean_tmp1 = mean(tri1.X(in.list_edges{i},:));
    mean_tmp2 = mean(tri2.X(in.list_edges{i},:));
    
    lapl_tri1 = tri1.X(i,:) - mean_tmp1;
    lapl_tri2 = tri2.X(i,:) - mean_tmp2;
    
    norm2_lapl_tri1 = sum( lapl_tri1.^2 );
    norm2_lapl_tri2 = sum( lapl_tri2.^2 );
    
    % TODO: bring out the tri1 - mean_tmp1
    F(i + in.n)     = in.sqrt_lambda * norm2_lapl_tri1;
    F(i + incr_row) = in.sqrt_lambda * norm2_lapl_tri2;
    
    jsi1 = (i - 1) * 3 + 1;
    jsi2 = jsi1 + incr_col;
    Js(i + in.n    , jsi1:jsi1+2) = (2 * in.sqrt_lambda).* lapl_tri1;%(in.sqrt_lambda * ( 1 / norm_lapl_tri1 ) ) .* lapl_tri1; 
    Js(i + incr_row, jsi2:jsi2+2) = (2 * in.sqrt_lambda).* lapl_tri2;%(in.sqrt_lambda * ( 1 / norm_lapl_tri2 ) ) .* lapl_tri2; 
    
    for j = 1:length(in.list_edges{i})
        edgein1 = (in.list_edges{i}(j)-1) * 3 + 1;
        edgein2 = edgein1 + incr_col;
        Js(i + in.n,     edgein1:edgein1+2) = Js(i + in.n,     jsi1:jsi1+2).*[in.scalar(i) in.scalar(i) in.scalar(i)];
        Js(i + incr_row, edgein2:edgein2+2) = Js(i + incr_row, jsi2:jsi2+2).*[in.scalar(i) in.scalar(i) in.scalar(i)];
    end

end

J = Js;
