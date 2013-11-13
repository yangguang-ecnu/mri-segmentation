function show_bound(BW,col) 

[B,L] = bwboundaries(BW,'noholes');
b = label2rgb(L, @jet, [1 0 0]);
if ~isempty(B)
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), col, 'LineWidth', 2)
    end 
end