magnitude = [0 2 3 4 5 6 7 8];
root = 'deform_';

for i=1:length(magnitude)
    disp(['magnitude ', num2str(magnitude(i))])
    save_name = strcat(root,num2str(magnitude(i)),'.mat');
    load_name = strcat('def_',num2str(magnitude(i)),'.mat');
    main_evaluation(save_name,load_name,dcmdir);
    
end