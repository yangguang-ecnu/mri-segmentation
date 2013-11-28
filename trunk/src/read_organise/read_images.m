function [images alpha_ch] = read_images(folder,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read and show the images from
%% a given folder.
%% Inputs: 1. folder -> name of folder containing files
%%         3. show ->  0/1 = show/not show the images
%%
%% Output: 1. images -> matrix containing all the imagesaquest
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dcm = dir(folder);

for i=3:size(dcm,1)
    
    if strcmp(dcm(i).name(end-7:end),'mask.png')
        
        filename = strcat(folder,dcm(i).name);
        [im map alpha] = imread(filename);

%         [r c] = find(im(:,:,1) == 0 & im(:,:,2) == 0 & im(:,:,3) == 0 );
%         
%         im(:,:,4) = zeros(size(im,1),size(im,2));
%         
%         for j=1:length(r)
%             im(r(j),c(j),1:3) = [255 255 255];
%             im(:,:,4) = 1;
%         end
%         images{i-2} = im(:,:,1:3).*[im(:,:,4);im(:,:,4);im(:,:,4)];

        images{i-2} = im;
        alpha_ch{i-2} = alpha;
        
        if show == 1
            h = imshow(images{i-2});title(['Image: ',dcm(i).name]);
            set(h,'AlphaData',alpha_ch{i-2});
            pause;
            close;
        end
    end
end

if show == 1
    close all;
end