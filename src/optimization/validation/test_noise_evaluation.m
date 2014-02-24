noise_mag = [0 1 2 3 4 5 6 7 8 9 10 11 15];
root = 'noise_';

load_name = 'def_4.mat';
load(load_name);

for i=1:length(noise_mag)
    
    disp(['------ Noise ', num2str(noise_mag(i)),'------------------------------------------------'])
    save_name = strcat(root,num2str(noise_mag(i)),'.mat');
    %% Apply noise
    opt_im_ax  = add_Rician_vol(opt_im_ax,  noise_mag(i));
    opt_im_sag = add_Rician_vol(opt_im_sag, noise_mag(i));
    opt_im_cor = add_Rician_vol(opt_im_cor, noise_mag(i));
    
    main_evaluation(save_name, load_name, dcmdir, opt_im_ax, opt_im_sag, opt_im_cor);
    
end