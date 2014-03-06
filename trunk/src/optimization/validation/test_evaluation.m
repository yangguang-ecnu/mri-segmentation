%% Change the folders in order to evaluate the different deformations %%
%% Here the deformation has been already computed %%

magnitude = [3 5 6 7 8 10 11 13 15];
root_save = 'deform_noise_';
root_load = 'deformations/def_';

for i=1:length(magnitude)
    
    disp('------------------------------------------------------------------------------------------------------')
    disp(['------------ Deformation Magnitude ', num2str(magnitude(i)),' -------------------------------------------'])
    disp('------------------------------------------------------------------------------------------------------')
    
    save_name = strcat(root_save,num2str(magnitude(i)),'.mat');
    load_name = strcat(root_load, num2str(magnitude(i)),'.mat'); %
    def_mag = load(load_name);
%     
%     im_ax  = add_Rician_vol(def_mag.opt_im_ax,  magnitude(i));
%     im_sag = add_Rician_vol(def_mag.opt_im_sag, magnitude(i));
%     im_cor = add_Rician_vol(def_mag.opt_im_cor, magnitude(i));

    main_evaluation(save_name, dcmdir, def_mag.opt_im_ax, def_mag.opt_im_sag, def_mag.opt_im_cor);
    
end