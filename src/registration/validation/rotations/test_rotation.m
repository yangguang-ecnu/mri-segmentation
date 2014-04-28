%% Rotation evaluation %%

fold = dir('rotations/');

for i=3:size(fold,1)
    
    filename = strcat(folder,fold(i).name);
    
    out = load(filename);
    disp(['-------- Rotation ',num2str(out.theta1),' ',num2str(out.theta2),' ',num2str(out.theta3), '-------------------------'])
    
    save_name = strcat('rot_',num2str(i-2),'.mat');
    
    main_evaluation(save_name, dcmdir, out.opt_im_ax, out.opt_im_sag, out.opt_im_cor);

end