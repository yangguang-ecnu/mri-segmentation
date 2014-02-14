% Creating vectors for x and y axes 
x = xmin:xstp:xmax; 
y = ymin:ystp:ymax; 
x=x/1000.;% renormalizes m in km, only for plotting 
y=y/1000.;% renormalizes m in km, only for plotting 
 
% Visualizing results 
 % Ploting RO as colormap 
 figure; 
% TK=ROT90(T0); 
TK = T0;
 surf(x,y,TK); 
 view (0,90); 
 
% Remember that you can use other plotting schemes, such as: 
 
light; lighting phong; view (50,30); shading interp; 
 
 caxis([tmin tmax]); 
 colorbar; 
 title(['Temperature, K For time step ',num2str(cycle)]) 
 xlabel('X, km') 
 ylabel('Y, km') 
