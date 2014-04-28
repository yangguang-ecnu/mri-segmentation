function [target_tri_ax, target_tri_sag, target_tri_cor] = initialize_optimization(  )

global source_tri_v

tic
%% initialization
mesh0 = [source_tri_v.X(:,1:2);source_tri_v.X(:,2:3);source_tri_v.X(:,1) source_tri_v.X(:,3)]; % [source_tri_v.X(:,1:2);source_tri_v.X(:,2:3)]

options = optimset('Display','iter','MaxIter', 100, 'TolFun', .1, 'PlotFcns',@optimplotfval);
[xfinal_tmp fval exitflag output] = fminunc(@myfun_unc_ortho_eval, mesh0, options);
 
xfinal(1:size(source_tri_v.X,1),:) = [xfinal_tmp(1:size(source_tri_v.X,1),:) source_tri_v.X(:,3)];
xfinal(1+size(source_tri_v.X,1):2*size(source_tri_v.X,1),:)   = [source_tri_v.X(:,1) xfinal_tmp(1+size(source_tri_v.X,1):2*size(source_tri_v.X,1),:)];
xfinal(1+2*size(source_tri_v.X,1):3*size(source_tri_v.X,1),:) = [ xfinal_tmp(1+2*size(source_tri_v.X,1):3*size(source_tri_v.X,1),1) source_tri_v.X(:,2) xfinal_tmp(1+2*size(source_tri_v.X,1):3*size(source_tri_v.X,1),2)];


target_tri_ax  = TriRep(source_tri_v.Triangulation, xfinal(1:size(source_tri_v.X,1),:));
target_tri_sag = TriRep(source_tri_v.Triangulation, xfinal(size(source_tri_v.X,1)+1:2*size(source_tri_v.X,1),:));
target_tri_cor = TriRep(source_tri_v.Triangulation, xfinal(2*size(source_tri_v.X,1)+1:end,:));

tim = toc