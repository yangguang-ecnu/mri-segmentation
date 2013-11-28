function [im_eq,h_eq,s,h1]=histogram_eq(I,h,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Histogram equalization
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(nargchk(1, 3, nargin));
if nargin == 2
    n = 256;
end
if isa(I,'double')
    A=1;
elseif isa(I,'uint8')
    A=255;
end
[rows cols] = size(I);
%h = histdouble(I);
s = zeros(1,n);
for k=1:n
    s(k) = sum(h(1:k));%./(rows*cols);   
end

%h1 = (s - min(s))./(rows*cols - min(s));

im_eq = (s(floor(I.*((n-1)/A)+1.5))-min(s))./(rows*cols-min(s));

h_eq = histdouble(im_eq);

h1 = zeros(1,n);
for k=1:n
    h1(k) = sum(h_eq(1:k));%./(rows*cols);   
end
  
