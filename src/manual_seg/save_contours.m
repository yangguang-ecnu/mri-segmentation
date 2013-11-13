
folder = 'seg2/image';

bladder = load('bladder/bladder_axial.mat');
bones = load('bones2/bones_axial.mat');

for i=1:size(views.axial,3)

    
    imshow(views.axial(:,:,i),[])
    hold on
    
    show_bound(bladder.I_texture(:,:,i),'r');
    show_bound(bones.I_texture(:,:,i),'g');
    
    if i < 10
        filename1 = strcat(folder,'00');
    elseif i >= 10
        filename1 = strcat(folder,'0');
    end
    
    filename2 = strcat(filename1,[num2str(i),'.png']);
    
    saveas(gcf,filename2,'png');
end