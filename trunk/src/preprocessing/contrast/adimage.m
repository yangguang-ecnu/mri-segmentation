function out_image = adimage(I,range_in, range_out,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%    ADJUST IMAGE
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%
disp(nargchk(1, 4, nargin));
if nargin == 1
    gamma = 1;
    range_in = [min(min(I)) max(max(I))];
    range_out = [0 1];
elseif nargin == 3
    gamma = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rows, cols, ch] = size(I);
% Check that the ranges are not empty
if isempty(range_in)
    range_in = [0 1];
end
if isempty(range_out)
    range_out = [0 1];
end

low_in = range_in(1);
high_in = range_in(2);
low_out = range_out(1);
high_out = range_out(2);

%%%%% Transform %%%%%
%% All the c_{i} are auxiliar matrix %%
c1 = (I <= low_in);
c2 = (I >= high_in);
c3 = c1|c2;
c4 = ( low_in < I < high_in);
out = (I-low_in).*((high_out-low_out)/(high_in-low_in)) + low_out;
c5 = c4.*out;
c6 = ~c3;
out_image = c6.*out + c1.*low_out + c2.*high_out;
out_image = out_image.^gamma;
figure;
subplot(1,2,1);imshow(I);title('Original image');
subplot(1,2,2);imshow(out_image);title('Adjust image');