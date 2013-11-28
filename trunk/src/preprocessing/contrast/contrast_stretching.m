function out=contrast_stretching(I,m,E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 1/(1 + (m/r)^{E})
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(nargchk(1, 3, nargin));
if nargin == 1
    m = (max(max(I))-min(min(I)))/2;
    E = 20;
end

out = 1./(1 + (m./I).^E);

figure;
subplot(1,2,1);imshow(I);title('Original image');
subplot(1,2,2);imshow(out);title('Contrast&Stretching image');