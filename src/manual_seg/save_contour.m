folder = 'bones2/image';


for i=1:size(views.axial,3)

    
    imshow(views.axial(:,:,i),[])
    hold on
    
    show_bound(I_texture(:,:,i),'g');

    if i < 10
        filename1 = strcat(folder,'00');
    elseif i >= 10
        filename1 = strcat(folder,'0');
    end
    
    filename2 = strcat(filename1,[num2str(i),'.png']);
    
    saveas(gcf,filename2,'png');
end