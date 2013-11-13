%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Contest create contours
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all
clc

[rows cols] = size(views.axial(:,:,1));

for k = 7:16%1:size(views.axial,3)
    
    figure;
    imshow(views.axial(:,:,k),[]);hold on
    
    gray(:,:,1) = im2double(views.axial(:,:,k));
    gray(:,:,2) = im2double(views.axial(:,:,k));
    gray(:,:,3) = im2double(views.axial(:,:,k));
    
    n = str2num(input('How many components ?  ','s'));

    %I_texture(:,:,k) = zeros(rows,cols);

    if n
        for i=n%1:n

            r = 0;

            while r==0
               r = str2num(input('Press 1 whenever you are ready  ','s'));
            end

             h = impoly;
            h_position = getPosition(h);
            [Vertices Lines] = contour(h_position);

            vert{k,i} = Vertices;
            lines{k,i} = Lines;

            tmp = poly2mask(vert{k,i}(:,1),vert{k,i}(:,2),rows,cols);

            tmp = imdilate(tmp,ones(3,3));

            I_texture(:,:,k) = I_texture(:,:,k) | tmp;


        end
    end
    
end