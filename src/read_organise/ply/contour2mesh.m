function [faces,vertices] = contour2mesh(mat_file,patient,ply_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  From a given set of contours (all the frames from one view) in a mat file
%%  from a given patient, compute the surface triangulation and save it in
%%  a ply file
%%
%%  NOTE: It is important to note that in the mat_file the contours are in 
%%        frame coordinates, rather than in the output ply file, the contour
%%        points are given in the RCS !
%%
%%  Inputs: 1. mat_file -> (string) of the path of the .mat file that contains
%%                         the structure 'vertex', a cell array where each
%%                         cell contains the vertex of the contour in one slice
%%          2. patient ->  (cell) with the DICOM info for each slice
%%          3. ply_file(optional) -> (string) path of the ply file to save the 
%%                                   vertex and faces information
%%
%%  Output: 1. faces -> Nfx3 matrix, where each row is a face and the columns
%%                              correspond to the vertex of this face
%%          2. vertices -> Nv x 3 matrix with rows number of contour points, 
%%                    and columns (x,y,z) coordinates
%%
%% EXECUTE
%% main_read;
%% [faces,vertices] = contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_sag2.mat',...
%%                         views.sagittal_info,'lassalas_sag.ply');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_file = 1;
if nargin < 3
    save_file = 0;
end

a = load(mat_file);

vert = a.vertex;

%% number of slices
% N = length(vert);
% length(patient)

x = [];
y = [];
z = [];

for i = 1:length(patient)
    
    if ~isempty(vert{i})
        M = zeros(4,4);
        
        M(1:3,4) = patient{i}.ImagePositionPatient;
        M(4,4) = 1;
        M(1:3,1) = patient{i}.ImageOrientationPatient(1:3).*patient{i}.PixelSpacing(1);
        M(1:3,2) = patient{i}.ImageOrientationPatient(4:6).*patient{i}.PixelSpacing(2);
        
        
        [Vertices Lines] = contour(vert{i});
        
        for k=1:size(Vertices,1)
            
            p = M*[Vertices(k,1) Vertices(k,2) 0 1]';
            
            x = [x p(1)];
            y = [y p(2)];
            z = [z p(3)];
            
        end
    end
    
end

var(:,1) = x(:)';
var(:,2) = y(:)';
var(:,3) = z(:)';

% dt = DelaunayTri(var);
% tetramesh(dt, 'FaceColor', 'cyan');

[t] = MyCrust(var);

faces = t;
vertices = var;

%% plot of the output triangulation
figure(2)
hold on
title('Output Triangulation','fontsize',14)
axis equal
trisurf(t,var(:,1),var(:,2),var(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
hold on

if save_file
    
    fileId = fopen(ply_file, 'wt');
    
    % writing the header
    fprintf(fileId, 'ply\n');
    fprintf(fileId, 'format ascii 1.0\n');
    fprintf(fileId, 'comment created by safir\n');
    
    % writing the elements
    fprintf(fileId, 'element vertex %d\n', length(x));
    % [x y z] x,y image coordinates; z the current slice
    fprintf(fileId, 'property float32 x\n');
    fprintf(fileId, 'property float32 y\n');
    fprintf(fileId, 'property float32 z\n');
    
    fprintf(fileId, 'element face %d\n', size(t,1));
    fprintf(fileId, 'property list int32 int32 vertex_indices\n');
    
    
    fprintf(fileId, 'end_header\n');
    
    % writing the vertexes
    
    for i=1:length(x)
        fprintf(fileId, '%f %f %f \n', x(i), y(i),z(i));
    end
    
    % writing the faces
    for i = 1 : size(t,1)
        fprintf(fileId, '3 %d %d %d \n', int32(t(i,1)), int32(t(i,2)), int32(t(i,3)));
    end
    fclose(fileId);
    
end
