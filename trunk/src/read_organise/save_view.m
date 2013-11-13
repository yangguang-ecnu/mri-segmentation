function save_view(volume,folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Save the given view (volume) in the given folder path
%%  NOTE: the path should not contain '/' at the end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(volume,3)
    
    folder = strcat(folder,'/image');
    
    if i < 10
        filename1 = strcat(folder,'00');
    elseif i >= 10
        filename1 = strcat(folder,'0');
    end
    
    filename2 = strcat(filename1,[num2str(i),'.png']);
    imwrite(convert2u8(volume(:,:,i)),filename2);
    
end