function [I_eq h_eq cdf h1] = histogram_eq_image(I,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Check channels of I
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows cols ch] = size(I);

if ch == 1
    [I_eq h_eq cdf h1] = histogram_eq(I,h);
elseif ch == 3
    [I_eq h_eq cdf h1] = histogram_eq_rgb(I,h);
end