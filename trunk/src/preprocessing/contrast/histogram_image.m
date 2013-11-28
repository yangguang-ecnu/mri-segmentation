function [h x n]= histogram_image(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Check the image size
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows cols ch] = size(I);

if ch == 1
    [h x]= histdouble(I);
    n=1;
elseif ch == 3
    [h x]= hist_rgb(I);
    n=3;
end