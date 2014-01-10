function exportPcbToPly( pcb, file_path )
%EXPORTPCBTOPLY exporting the simulated pcb to a ply file
%   Detailed explanation goes here

fileId = fopen(file_path, 'wt');

% writing the header
fprintf(fileId, 'ply\n');
fprintf(fileId, 'format ascii 1.0\n');
fprintf(fileId, 'comment created by safir\n');
fprintf(fileId, 'comment id look up table\n');
for i = 1 : length(pcb.idList)
   fprintf(fileId, 'comment id %d value %s \n', i, pcb.idList{i});
end

fprintf(fileId, 'comment part look up table\n');
for i = 1 : length(pcb.partList)
   fprintf(fileId, 'comment part %d value %s \n', i, pcb.partList{i});
end

fprintf(fileId, 'comment texture look up table\n');
for i = 1 : length(pcb.textureList)
   fprintf(fileId, 'comment texture %d value %s \n', i, pcb.textureList{i});
end

% writing the elements
[nvertex ~] = size(pcb.vertex);
fprintf(fileId, 'element vertex %d\n', nvertex);
% [x y z] position
fprintf(fileId, 'property float32 x\n');
fprintf(fileId, 'property float32 y\n');
fprintf(fileId, 'property float32 z\n');

[nfaces ~] = size(pcb.face);
fprintf(fileId, 'element face %d\n', nfaces);
fprintf(fileId, 'property list int32 int32 vertex_indices\n');
fprintf(fileId, 'property uchar red\n');
fprintf(fileId, 'property uchar green\n');
fprintf(fileId, 'property uchar blue\n');
fprintf(fileId, 'property int label\n');

fprintf(fileId, 'end_header\n');

% writing the vertexes
for i = 1 : nvertex
   fprintf(fileId, '%f %f %f \n', pcb.vertex(i, 1), pcb.vertex(i, 2), pcb.vertex(i, 3));
end

% writing the faces
for i = 1 : nfaces
   fprintf(fileId, '3 %d %d %d %u %u %u %02d%02d%02d \n', pcb.face(i, 1) - 1, pcb.face(i, 2) - 1, pcb.face(i, 3) - 1, pcb.rgb(i, 1), pcb.rgb(i, 2), pcb.rgb(i, 3), pcb.id(i), pcb.part(i), pcb.texture(i));
end

fclose(fileId);

end