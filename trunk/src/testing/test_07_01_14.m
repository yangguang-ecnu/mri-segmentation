clc;
close all;

%% Set up the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('set_up','var')
    set_parameters;
end

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------- Start Optimization  -----')

size1 = size(var_cell1,3);
size2 = size(var_cell1,2);
size3 = size(var_cell1,1);

global sub_1
global sub_2
global sub_3
global current_tr
global c2b_coord

[sub_1, sub_2, sub_3] = ind2sub([size1 size2 size3],1:size(var_array1,1));
current_tr = tsearchn(source_tri.X,source_tri.Triangulation,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
c2b_coord = cartToBary(source_tri,current_tr,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % barycentric coordinates of the points wrt source tri


lbv = [-20 -20 -20];
ubv = [20 20 20];

lb = repmat(lbv,2*size(source_tri.X,1),1) + [source_tri.X;source_tri.X];
ub = repmat(ubv,2*size(source_tri.X,1),1) + [source_tri.X;source_tri.X];
% lb = [];
% ub = [];

opts = optimset('Jacobian','on','Display','iter','MaxIter', 20);

tic
[in out] = preparing;
tim = toc

tic
% [xfinal fval exitflag output] = lsqnonlin(@(t)myfun4(t,in),[source_tri.X;source_tri.X],lb,ub,opts);
[xfinal fval exitflag output] = lsqnonlin(@(t)myfun5(t,in,out),[source_tri.X;source_tri.X],lb,ub,optimset('Display','iter','MaxIter', 20));

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




