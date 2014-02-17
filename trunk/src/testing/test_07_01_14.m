clc;
close all;

%% Set up the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('set_up','var')
    set_parameters;
end

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------- Start Optimization  -----')

global source_tri
global var_array1
global var_cell1
global var_array2
global var_cell2
global var_array3
global var_cell3

disp('--------- Preparing first intersection: axial/sagittal  -----')
size1_int1 = size(var_cell1,3);
size2_int1 = size(var_cell1,2);
size3_int1 = size(var_cell1,1);

global sub_1_int1
global sub_2_int1
global sub_3_int1
global current_tr_int1
global c2b_coord_int1


[sub_1_int1, sub_2_int1, sub_3_int1] = ind2sub([size1_int1 size2_int1 size3_int1],1:size(var_array1,1));
current_tr_int1 = tsearchn(source_tri.X,source_tri.Triangulation,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int1  = cartToBary(source_tri,current_tr_int1,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % barycentric coordinates of the points wrt source tri

disp('--------- Preparing second intersection: axial/coronal  -----')
size1_int2 = size(var_cell2,3);
size2_int2 = size(var_cell2,2);
size3_int2 = size(var_cell2,1);

global sub_1_int2
global sub_2_int2
global sub_3_int2
global current_tr_int2
global c2b_coord_int2


[sub_1_int2, sub_2_int2, sub_3_int2] = ind2sub([size1_int2 size2_int2 size3_int2],1:size(var_array2,1));
current_tr_int2 = tsearchn(source_tri.X,source_tri.Triangulation,[var_array2(:,1) var_array2(:,2) var_array2(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int2  = cartToBary(source_tri,current_tr_int2,[var_array2(:,1) var_array2(:,2) var_array2(:,3)]); % barycentric coordinates of the points wrt source tri


disp('--------- Preparing third intersection: sagittal/coronal  -----')
size1_int3 = size(var_cell3,3);
size2_int3 = size(var_cell3,2);
size3_int3 = size(var_cell3,1);

global sub_1_int3
global sub_2_int3
global sub_3_int3
global current_tr_int3
global c2b_coord_int3


[sub_1_int3, sub_2_int3, sub_3_int3] = ind2sub([size1_int3 size2_int3 size3_int3],1:size(var_array3,1));
current_tr_int3 = tsearchn(source_tri.X,source_tri.Triangulation,[var_array3(:,1) var_array3(:,2) var_array3(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord_int3  = cartToBary(source_tri,current_tr_int3,[var_array3(:,1) var_array3(:,2) var_array3(:,3)]); % barycentric coordinates of the points wrt source tri


lbv = [-20 -20 -20];
ubv = [ 20  20  20];

lb = repmat(lbv,2*size(source_tri.X,1),1) + [source_tri.X;source_tri.X];
ub = repmat(ubv,2*size(source_tri.X,1),1) + [source_tri.X;source_tri.X];
% lb = [];
% ub = [];

opts = optimset('Jacobian','on','Display','iter','MaxIter', 20);

tic
% [in out] = preparing;
preparing_eval;
% [in out] = preparing4;
tim = toc

tic
%% initialization
mesh0 = [source_tri.X(:,1:2);source_tri.X(:,2:3);source_tri.X(:,1) source_tri.X(:,3)]; % [source_tri.X(:,1:2);source_tri.X(:,2:3)]

% opt4 = optimset('Algorithm','levenberg-marquardt','Jacobian','on','DerivativeCheck','on','GradConstr','on');
% [xfinal fval exitflag output] = lsqnonlin(@(t)myfun4(t,in),[source_tri.X;source_tri.X],lb,ub,opt4);
% [xfinal fval exitflag output] = lsqnonlin(@(t)myfun5(t,in,out),[source_tri.X;source_tri.X],lb,ub,optimset('Display','iter','MaxIter', 40));
% [xfinal fval exitflag output] = fminunc(@myfun_unc,[source_tri.X;source_tri.X],optimset('Display','iter','MaxIter', 40));
% [xfinal_tmp fval exitflag output] = fminunc(@myfun_unc_ortho, mesh0, optimset('Display','iter','MaxIter', 40));

% options = optimset('LargeScale','off','Jacobian','on','Display','iter','DerivativeCheck','on');
options = optimset('Display','iter','MaxIter', 40);
[xfinal_tmp fval exitflag output] = fminunc(@myfun_unc_ortho, mesh0,options);
% 
xfinal(1:size(source_tri.X,1),:) = [xfinal_tmp(1:size(source_tri.X,1),:) source_tri.X(:,3)];
xfinal(1+size(source_tri.X,1):2*size(source_tri.X,1),:)   = [source_tri.X(:,1) xfinal_tmp(1+size(source_tri.X,1):2*size(source_tri.X,1),:)];
xfinal(1+2*size(source_tri.X,1):3*size(source_tri.X,1),:) = [ xfinal_tmp(1+2*size(source_tri.X,1):3*size(source_tri.X,1),1) source_tri.X(:,2) xfinal_tmp(1+2*size(source_tri.X,1):3*size(source_tri.X,1),2)];

%% Define the initial mesh as a vector
% mesh_vector_tmp = [source_tri.X;source_tri.X]';
% mesh_vector = mesh_vector_tmp(:);
% 
% [xfinal_v, fval] = optimizationWithDE(1, 6*size(source_tri.X,1),[],[],[], [], [], [], [], []);
% xfinal = reshape(xfinal_v',2*size(source_tri.X,1),3) +  [source_tri.X;source_tri.X];

tim = toc


disp('--------- Plot the data points  -----')

%Plot the data points  

plot3(var_array1(:,1),var_array1(:,2),var_array1(:,3),'r*','MarkerSize',2);hold on
plot3(var_array2(:,1),var_array2(:,2),var_array2(:,3),'b*','MarkerSize',2);hold on
plot3(var_array3(:,1),var_array3(:,2),var_array3(:,3),'k*','MarkerSize',2);hold on


% Plot the points
% figure(1);
% for i = 1:length(n_points)
%     plot3(sol(i,1),sol(i,2),sol(i,3),'r*');hold on
%     plot3(var_array(n_points(i),1),var_array(n_points(i),2),var_array(n_points(i),3),'g*');hold on
%     line([sol(i,1),var_array(n_points(i),1)],[sol(i,2),var_array(n_points(i),2)],[sol(i,3),var_array(n_points(i),3)]);hold on
% end


%% Clear memory

varlist = {'X_ax','Y_ax','Z_ax','X_sag','Y_sag','Z_sag','M','rows','cols','total_ax','total_sag','N1','N2','ax','sag','x','y','z','options','i','j','k','A1','b1','Aeq1','beq1','A2','b2','Aeq2','beq2',...
           'A','b','Aeq','beq','x0','lambda','V1','nr','ne','vd','ind_tmp','ind_tmp2','ind_tmp3','ind_tmp4','ind_tmp5','ind_tmp6','ind_tmp7','ind_tmp8','tmp_v1','tmp_v',...
           'new_i','new_j','new_k','neig','A_c','b_c','Aeq_c','beq_c','bb','nx','ny','nz','tmp','s2ind','l_x','l_y','l_z','vari','dt','trep'};
clear(varlist{:})
clear varlist




