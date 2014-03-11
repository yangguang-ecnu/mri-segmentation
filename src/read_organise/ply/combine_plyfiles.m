function [faces,vertices] = combine_plyfiles(ply_ax,ply_sag,ply_cor,combine_ply, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Given the 3 ply files of the contours of the 3 different views 
%%  combine all the vertex and calculate the new mesh.
%%
%%  NOTE: the coordinate points in the ply files (x,y,z) are in the 3D space
%%        i.e, in the RCS !
%%
%%  Inputs:  1. ply_ax -> (string) ply file containing the vertex 
%%                        of the segmented organ, axial view 
%%           2. ply_sag -> (string) ply file containing the vertex 
%%                        of the segmented organ, sagittal view 
%%           3. ply_cor -> (string) ply file containing the vertex 
%%                        of the segmented organ, coronal view 
%%           4. combine_ply(optional) -> (string) path of the ply file to save the 
%%                                   vertex and faces information
%%           5. show    -> if 1 shows the calculate surface, 0 (default) otherwise  
%% Outputs:  1. faces -> Nfx3 matrix, where each row is a face and the columns
%%                              correspond to the vertex of this face
%%           2. vertices -> Nv x 3 matrix with rows number of contour points, 
%%                    and columns (x,y,z) coordinates
%%
%% EXECUTE
%%
%% contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_cor2.mat',views.coronal_info,'lassalas_cor.ply');
%% contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_ax2.mat',views.axial_info,'lassalas_ax.ply');
%% contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_sag2.mat',views.sagittal_info,'lassalas_sag.ply');
%% combine_plyfiles('lassalas_ax.ply','lassalas_sag.ply','lassalas_cor.ply');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_file = 1;

%% Check the number of inputs 
if nargin == 5
    [va,fa] = read_ply(ply_ax);
    [vs,fs] = read_ply(ply_sag);
    [vc,fc] = read_ply(ply_cor);
    
    var(:,1) = [va(:,1);vs(:,1);vc(:,1)];
    var(:,2) = [va(:,2);vs(:,2);vc(:,2)];
    var(:,3) = [va(:,3);vs(:,3);vc(:,3)];
end

if nargin == 4
    [va,fa] = read_ply(ply_ax);
    [vs,fs] = read_ply(ply_sag);
    [vc,fc] = read_ply(ply_cor);
    
    var(:,1) = [va(:,1);vs(:,1);vc(:,1)];
    var(:,2) = [va(:,2);vs(:,2);vc(:,2)];
    var(:,3) = [va(:,3);vs(:,3);vc(:,3)];
    show = 0;
end

if nargin == 3 
    save_file = 0;
    show = 0;
    [va,fa] = read_ply(ply_ax);
    [vs,fs] = read_ply(ply_sag);
    [vc,fc] = read_ply(ply_cor);
    
    var(:,1) = [va(:,1);vs(:,1);vc(:,1)];
    var(:,2) = [va(:,2);vs(:,2);vc(:,2)];
    var(:,3) = [va(:,3);vs(:,3);vc(:,3)];
end

if nargin == 2
    save_file = 0;
    show = 0;
    [va,fa] = read_ply(ply_ax);
    [vs,fs] = read_ply(ply_sag);
    
    var(:,1) = [va(:,1);vs(:,1)];
    var(:,2) = [va(:,2);vs(:,2)];
    var(:,3) = [va(:,3);vs(:,3)];
end

if nargin == 1
    save_file = 0;
    show = 0;
    [va,fa] = read_ply(ply_ax);
    
    var(:,1) = va(:,1);
    var(:,2) = va(:,2);
    var(:,3) = va(:,3);
end


[t] = MyCrust(var);

FV.vertices = var;
FV.faces = t;

%% Plot of the output triangulation
if show
    figure
    hold on
    title('Output Triangulation','fontsize',14)
    axis equal
    subplot(121);trisurf(t,var(:,1),var(:,2),var(:,3),'facecolor','c','edgecolor','b');camlight
    hold on
end

%% Plot the smooth triangulation

FV2 = smoothpatch(FV,1,5);

faces = FV2.faces;
vertices = FV2.vertices;

if show
    title('Output Smooth Triangulation','fontsize',14)
    axis equal
    subplot(122);trisurf(FV2.faces,FV2.vertices(:,1),FV2.vertices(:,2),FV2.vertices(:,3),'facecolor','c','edgecolor','b','edgealpha',0,'facelighting','flat');camlight
    hold on
end

%% Save the contours in a ply file 
if save_file
    
    fileId = fopen(combine_ply, 'wt');
    
    % writing the header
    fprintf(fileId, 'ply\n');
    fprintf(fileId, 'format ascii 1.0\n');
    fprintf(fileId, 'comment created by safir\n');
    
    % writing the elements
    fprintf(fileId, 'element vertex %d\n', size(FV2.vertices,1));%var
    
    fprintf(fileId, 'property float32 x\n');
    fprintf(fileId, 'property float32 y\n');
    fprintf(fileId, 'property float32 z\n');
    
%     fprintf(fileId, 'element face %d\n', size(FV2.faces,1));%t
%     fprintf(fileId, 'property list int32 int32 vertex_indices\n');
%     
%     fprintf(fileId, 'property uchar red\n');
%     fprintf(fileId, 'property uchar green\n');
%     fprintf(fileId, 'property uchar blue\n');
    
    fprintf(fileId, 'end_header\n');
    
    % writing the vertexes
    
    for i=1:size(FV2.vertices,1)
        fprintf(fileId, '%f %f %f \n', FV2.vertices(i,1),FV2.vertices(i,2),FV2.vertices(i,3));
    end
    
%     % writing the faces
%     for i = 1 : size(FV2.faces,1)
%         fprintf(fileId, '3 %d %d %d %u %u %u\n', int32(FV2.faces(i,1)), int32(FV2.faces(i,2)), int32(FV2.faces(i,3)), 0, 255, 255);
%     end
    fclose(fileId);
    
end