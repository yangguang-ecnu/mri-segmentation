function preparing_eval



global vol_ax
global vol_sag
global vol_cor

global list_edges
global neig_ax
global neig_sag
global neig_cor


global scalar
for i = 1:length(list_edges)
    scalar(i) = -1/length(list_edges{i});
end


%% Neighborhood
size_neig = [5 5];

for k = 1:size(vol_ax,3);
    neig_ax{k}  = colfilt( vol_ax(:,:,k), size_neig,'sliding',@mean);
end

for k = 1:size(vol_sag,3);
    neig_sag{k} = colfilt( vol_sag(:,:,k),size_neig,'sliding',@mean);
end

for k = 1:size(vol_cor,3);
    neig_cor{k} = colfilt( vol_cor(:,:,k),size_neig,'sliding',@mean);
end



