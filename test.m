
tmp = views.axial;

for k = 1:size(views.axial,3)
    for j = 1:512
        a(k,j) = views.axial(250,j,k);
    end
end
a = imresize(a,[512,512],'bilinear');
%subplot(221);imshow(views.axial(:,:,1),[]);
figure;imshow(a,[]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% cla
% 
% 
% cor = 7:11;
% 
% 
% for ax = 3:4:21
%     
%     aa(:,:,1) = convert2u8(views.axial(:,:,ax));
%     aa(:,:,2) = convert2u8(views.axial(:,:,ax));%zeros(size(handles.vol_axial(:,:,ax)));
%     aa(:,:,3) = convert2u8(views.axial(:,:,ax));%zeros(size(handles.vol_axial(:,:,ax)));
% 
% 
% 
% 
% 
%     %% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     M = zeros(4,4);
% 
%     M(1:3,4) = views.axial_info{ax}.ImagePositionPatient;
%     M(4,4) = 1;
%     M(1:3,1) = views.axial_info{ax}.ImageOrientationPatient(1:3).*views.axial_info{ax}.PixelSpacing(1);
%     M(1:3,2) = views.axial_info{ax}.ImageOrientationPatient(4:6).*views.axial_info{ax}.PixelSpacing(2);
% 
%     % cross_f = cross(handles.views.axial_info{ax}.ImageOrientationPatient(4:6),handles.views.axial_info{ax}.ImageOrientationPatient(1:3))
%     % 
%     % M(1:3,3) = cross_f .* handles.views.axial_info{ax}.SpacingBetweenSlices;
% 
%     i = [0 511];
%     j = [0 511];
% 
%     x = [];
%     y = [];
%     z = [];
% 
% 
%     for k=1:length(i)
%         for l=1:length(j)
%         p = M*[j(k) i(l) double(ax) 1]';
% 
%             x = [x p(1)];
%             y = [y p(2)];
%             z = [z p(3)];
% 
%         end
%     end
% 
%     P1 = [x(1) y(1) z(1)];
%     P2 = [x(2) y(2) z(2)];
%     P3 = [x(3) y(3) z(3)];
%     P4 = [x(4) y(4) z(4)];
% 
%     x = [P1(1) P4(1) P2(1)];   %# [xorigin xA xB] coordinates in 3-D space
%     y = [P1(2) P4(3) P2(2)];   %# [yorigin yA yB] coordinates in 3-D space
%     z = [P1(3) P4(3) P2(3)];   %# [zorigin zA zB] coordinates in 3-D space
% 
%     origin = [512 1];  %# Vertex of triangle in image space
%     U = [0 511];       %# Vector from origin to point A in image space
%     V = [-511 511]; 
%     img = aa;  %# Sample image for texture map
% 
%     A = origin+U;  %# Point A
%     B = origin+V;  %# Point B
%     C = B-U;     %# Point C
% 
%     [nRows,nCols,nPages] = size(img);  %# Image dimensions
%     inputCorners = [origin; ...        %# Corner coordinates of input space
%                     A; ...
%                     B; ...
%                     C];
% 
%     outputCorners = [origin; ...      %# Corner coordinates of output space
%                      A; ...
%                      B; ...
%                      C];            
% 
%     tform = maketform('projective',...  %# Make the transformation structure
%                       inputCorners,...
%                       outputCorners);
%     triTexture = imtransform(img,tform,'bicubic',...  %# Transform the image
%                              'xdata',[1 nCols],...
%                              'ydata',[1 nRows],...
%                              'size',[nRows nCols]);
%     x = [P3(1) P4(1) P2(1) P1(1)];   %# [xorigin xA xB] coordinates in 3-D space
%     y = [P3(2) P4(2) P2(2) P1(2)];   %# [yorigin yA yB] coordinates in 3-D space
%     z = [P3(3) P4(3) P2(3) P1(3)];   %# [zorigin zA zB] coordinates in 3-D space
%     % index = [3 4; 2 1];  %# Index used to create 2-by-2 surface coordinates
%     %index = [1 4;2 3];
%     index = [4 1;3 2];
%     X1 = x(index);      %# x coordinates of surface
%     Y1 = y(index);       %# y coordinates of surface
%     Z1 = z(index);        %# z coordinates of surface
% 
% 
% 
%     hSurface = surf(X1,Y1,Z1,triTexture,...          %# Plot texture-mapped surface
%                     'FaceColor','texturemap',...
%                     'EdgeColor','none');hold on
%     axis equal 
% 
% end
% 
% %% Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for cor = 3:5:23
%     
%     cc(:,:,3) = convert2u8(views.coronal(:,:,cor));
%     cc(:,:,1) = convert2u8(views.coronal(:,:,cor));%zeros(size(handles.vol_coronal(:,:,cor)));
%     cc(:,:,2) = convert2u8(views.coronal(:,:,cor));%zeros(size(handles.vol_coronal(:,:,cor)));
% 
%     M = zeros(4,4);
% 
%     M(1:3,4) = views.coronal_info{cor}.ImagePositionPatient;
%     M(4,4) = 1;
%     M(1:3,1) = views.coronal_info{cor}.ImageOrientationPatient(1:3).*views.coronal_info{cor}.PixelSpacing(1);
%     M(1:3,2) = views.coronal_info{cor}.ImageOrientationPatient(4:6).*views.coronal_info{cor}.PixelSpacing(2);
% 
% 
% 
%     i = [0 511];
%     j = [0 511];
%     x = [];
%     y = [];
%     z = [];
%     for k=1:2
%         for l=1:2
%         p = M*[j(k) i(l) 0 1]';
% 
%             x = [x p(1)];
%             y = [y p(2)];
%             z = [z p(3)];
%         end
%     end
% 
% 
%     P1 = [x(1) y(1) z(1)];
%     P2 = [x(2) y(2) z(2)];
%     P3 = [x(3) y(3) z(3)];
%     P4 = [x(4) y(4) z(4)];
% 
%     x = [P1(1) P4(1) P2(1)];   %# [xorigin xA xB] coordinates in 3-D space
%     y = [P1(2) P4(3) P2(2)];   %# [yorigin yA yB] coordinates in 3-D space
%     z = [P1(3) P4(3) P2(3)];   %# [zorigin zA zB] coordinates in 3-D space
% 
%     origin = [512 1];  %# Vertex of triangle in image space
%     U = [0 511];       %# Vector from origin to point A in image space
%     V = [-511 511]; 
%     img = cc;  %# Sample image for texture map
% 
%     A = origin+U;  %# Point A
%     B = origin+V;  %# Point B
%     C = B-U;     %# Point C
% 
%     [nRows,nCols,nPages] = size(img);  %# Image dimensions
%     inputCorners = [origin; ...        %# Corner coordinates of input space
%                     A; ...
%                     B; ...
%                     C];
% 
%     outputCorners = [origin; ...      %# Corner coordinates of output space
%                      A; ...
%                      B; ...
%                      C];            
% 
%     tform = maketform('projective',...  %# Make the transformation structure
%                       inputCorners,...
%                       outputCorners);
%     triTexture = imtransform(img,tform,'bicubic',...  %# Transform the image
%                              'xdata',[1 nCols],...
%                              'ydata',[1 nRows],...
%                              'size',[nRows nCols]);
%     x = [P3(1) P4(1) P2(1) P1(1)];   %# [xorigin xA xB] coordinates in 3-D space
%     y = [P3(2) P4(2) P2(2) P1(2)];   %# [yorigin yA yB] coordinates in 3-D space
%     z = [P3(3) P4(3) P2(3) P1(3)];   %# [zorigin zA zB] coordinates in 3-D space
% 
%     index = [4 1;3 2];  %# Index used to create 2-by-2 surface coordinates
%     X1 = x(index);      %# x coordinates of surface
%     Y1 = y(index);       %# y coordinates of surface
%     Z1 = z(index);        %# z coordinates of surface
% 
% 
%     hSurface = surf(X1,Y1,Z1,triTexture,...          %# Plot texture-mapped surface
%                     'FaceColor','texturemap',...
%                     'EdgeColor','none');
%     axis equal
% 
% end
% 
% view(3)
% grid off
% axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [optimizer,metric] = imregconfig('multimodal');
% 
% tic
% movingRegisteredDefault = imregister(ax, ax2, 'affine', optimizer, metric);
% timeDefault = toc
% 
% figure, imshowpair(movingRegisteredDefault, ax2)
% title('A: Default registration')


% % % Compile the c-code
%   mex src/preprocessing/biasfield/Kroon/BCFCM3D.c -v
%   
% % Load test volume
%   D = views.axial;
%   
% % Convert to single
%   D = im2single(D);
% 
% % Class prototypes (means)
%   v = [0 0.12 0.18];
%   
% % Do the fuzzy clustering
%   [B,U] = BCFCM3D(D,v,struct('maxit',5,'epsilon',1e-5,'sigma',1));
%   
%   cl_vol = classes2gray(U);
% 
% % Show results
%   figure, 
%   subplot(2,2,1), imshow(D(:,:,15),[]), title('A slice of input volume');
%   subplot(2,2,2), imshow(squeeze(U(:,:,15,1:3))), title('3 classes of the partition matrix');
%   subplot(2,2,3), imshow(B(:,:,15),[]), title('Estimated biasfield');
%   subplot(2,2,4), imshow(D(:,:,15)-B(:,:,15),[]), title('Corrected slice');

% X = 0:65535;
% 
% uterus = normpdf(X,303.001,126.741);
% bladder = normpdf(X,638.03,289.75);
% vagina  = normpdf(X,888.18,428.166);
% rectum  = normpdf(X,114.589,161.649);
% 
% figure;
% plot(X,uterus,'g');hold on
% plot(X,bladder,'b');hold on
% plot(X,vagina,'r');hold on
% plot(X,rectum,'y');hold on
% axis([0 4000  0 .0035]);