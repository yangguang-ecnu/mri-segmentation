function image_mat = cellimages2mat(image_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = length(image_cell);

for i = 1:p
    if size(image_cell{i},3) ~= 1
        image_mat(:,:,i) = rgb2ind(image_cell{i},3);
    else
        image_mat(:,:,i) = image_cell{i};
    end
end