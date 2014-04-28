 % Compile the c-code
  mex src/preprocessing/biasfield/BCFCM2D.c -v;
 % Load test image
   %Y=im2double(imread('test_biasfield_noise.png'));
   Y = double(views.sagittal(:,:,25));
   level = max(Y(:))/100;
   M = 3;
   al = 1;
   %Y = ornlm(Y,M+2,al,level);
 % Class prototypes (means)
   v = [10;30;50];
 % Do the fuzzy clustering
   [B,U]=BCFCM2D(Y,v,struct('maxit',15,'epsilon',1e-5));
 % Show results
   figure, 
   subplot(2,2,1), imshow(Y,[]), title('input image');
   subplot(2,2,2), imshow(U,[]), title('Partition matrix');
   subplot(2,2,3), imshow(B,[]), title('Estimated biasfield');
   subplot(2,2,4), imshow(Y-B,[]), title('Corrected image');
