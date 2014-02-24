%% Compute graph between error and magnitude ....

mag = [0 2 3 5 6 7 8];

x = [];
y = [];

for i=1:length(mag)
    
    er = load(strcat('error_deform_',num2str(mag(i)),'.mat'));
    x = [x mag(i)];
    y = [y er.error];
    
end

figure;plot(x, y, 'ko');axis([-1 max(x)+1 0 max(y)+2]);