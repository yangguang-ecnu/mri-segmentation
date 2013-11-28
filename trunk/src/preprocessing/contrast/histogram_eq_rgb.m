function [I_eq h_eq cdf h1] = histogram_eq_rgb(I,h,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I_eq(:,:,1) h_eq(:,1) cdf(:,1) h1(:,1)] = histogram_eq(I(:,:,1),h(:,1));
[I_eq(:,:,2) h_eq(:,2) cdf(:,2) h1(:,2)] = histogram_eq(I(:,:,2),h(:,2));
[I_eq(:,:,3) h_eq(:,3) cdf(:,3) h1(:,3)] = histogram_eq(I(:,:,3),h(:,3));