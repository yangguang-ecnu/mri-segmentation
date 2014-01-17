%% Naive example 
%% Define a cube and compute the triangulation 

% cube(:,:,1) = [0 0 1;0 1 1;0 2 1; ...
%                1 0 1;1 1 1;1 2 1;...
%                2 0 1;2 1 1;2 2 1]; ...
% cube(:,:,2) = [0 0 0;0 1 0;0 2 0; ...
%                1 0 0;1 1 0;1 2 0; ...
%                2 0 0;2 1 0;2 2 0];
figure          
cube{1,1,1} = [0 0 1];cube{1,1,2} = [0 1 1];cube{1,1,3} = [0 2 1];
cube{1,2,1} = [1 0 1];cube{1,2,2} = [1 1 1];cube{1,2,3} = [1 2 1];
cube{1,3,1} = [2 0 1];cube{1,3,2} = [2 1 1];cube{1,3,3} = [2 2 1];

cube{2,1,1} = [0 0 0];cube{2,1,2} = [0 1 0];cube{2,1,3} = [0 2 0];
cube{2,2,1} = [1 0 0];cube{2,2,2} = [1 1 0];cube{2,2,3} = [1 2 0];
cube{2,3,1} = [2 0 0];cube{2,3,2} = [2 1 0];cube{2,3,3} = [2 2 0];

ax = 2;
sag = 3;
t = 3;

FV2.vertices = [];
FV2.faces = [];
tetra = [];
for i=1:ax-1%size(cube,3)
    for j=1:sag-1%size(cube(:,:,i),1)
        
        %vari = [];
        %vari = [vari;cube{i,j,1}];
        
        for k=1:t-1
            
            s2ind =  k + t*(j-1 + sag*(i-1));
            vert_mat(i,j,k) = s2ind;%ijk2index(i,j,k,ax,sag);
            
            vari = [cube{i,j,k};cube{i,j,k+1};cube{i+1,j,k};cube{i+1,j+1,k};cube{i,j+1,k};cube{i+1,j+1,k+1};cube{i,j+1,k+1};cube{i+1,j,k+1}];
            
%             FV2.vertices = [FV2.vertices;vari];
%             
%             dt =        [vert_mat(i,j,k) vert_mat(i,j+1,k)   vert_mat(i+1,j+1,k) 1;...
%                          vert_mat(i,j,k) vert_mat(i+1,j+1,k) vert_mat(i+1,j,k) 1; ...
%                          vert_mat(i,j,k) vert_mat(i,j,k+1) vert_mat(i,j+1,k+1) 1; ... 
%                          vert_mat(i,j,k) vert_mat(i,j+1,k) vert_mat(i,j+1,k+1) 1; ...
%                          vert_mat(i,j,k) vert_mat(i,j,k+1) vert_mat(i,j+1,k) 1; ...
%                          vert_mat(i,j,k) vert_mat(i+1,j,k) vert_mat(i,j+1,k) 1; ...
%                          vert_mat(i,j,k) vert_mat(i+1,j+1,k+1) vert_mat(i+1,j,k) 1];
%             FV2.faces = [FV2.faces; dt];
            %dt = DelaunayTri(vari);%
            tmp = calculate_tetrahedron(i,j,k,t,sag);
            
            tetra = [tetra;tmp];
%             dt.Triangulation
             %FV2.faces = dt.Triangulation;
             %FV2.vertices = vari;
            %trisurf(dt,vari(:,1),vari(:,2),vari(:,3),'facecolor','r','facealpha',.3,'edgecolor','b','facelighting','flat');hold on;%camlight
            %trisurf(FV2.faces,FV2.vertices(:,1),FV2.vertices(:,2),FV2.vertices(:,3),'facecolor','r','facealpha',.3,'edgecolor','b','facelighting','flat');hold on;%camlight
            plot3(cube{i,j,k}(1), cube{i,j,k}(2), cube{i,j,k}(3),'g*');hold on
            
            tetra_d.X = vari;
            tetra_d.Triangulation = [1,3,8,4;2,1,8,4;2,5,1,4;2,6,5,4;2,8,6,4;2,7,5,6];
            
            %tetramesh(tmp);hold on
            
            
            
        end
                
    end
end





