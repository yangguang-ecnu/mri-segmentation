function histogram(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the histogram of the image I
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = adimage(I,[],[0 255]);
[rows cols ch] = size(J);
h = zeros(1,256);
for i=1:rows
    for j=1:cols
        h(floor(J(i,j))+1) = h(floor(J(i,j))+1) + 1;
    end
end

x = [0:3:255];
bar(h,'b');title('HISTOGRAM');
%x = zeros(1,256);
%x(1:5:256) = h(1:5:256);
%stem(x,h(1:5:256),'fill');title('HISTOGRAM');
axis([0 255 0 max(max(h))])



colormap(gray(256))

colorbar('Location','SouthOutside','XTickLabel',{  });


