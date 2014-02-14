maxcycle = 10; % Defining the number of cycles for a simulation

% Defining grid in the horizontal direction
xmin=0.0; xmax=3000.0; xnum=99;
xstp=(xmax-xmin)/(xnum-1);
ymin=0.0; ymax=3000.0; ynum=99;
ystp=(ymax-ymin)/(ynum-1);

% Defining Material Properties
k=3.0; % Thermal conductivity, W/m/K
ro=3000.0; %Density, kg/m^3
cp=1000.0; % Isobaric heat capacity, J/kg
kappa=k/ro/cp; % Thermal diffusivity, m^2/s

% Defining Initial Temperature structure
tmin=500;
tmax=1000;
for xn = 1:xnum
    for yn = 1:ynum
        T0(xn,yn)=tmin; % Background Temperature, K
        if ((xn-1)/(xnum-1)>0.60 && (xn-1)/(xnum-1)<0.80 && (yn-1)/(ynum-1)>0.60 && (yn-1)/(ynum-1)<0.80)
            T0(xn,yn)=tmax; % Rectangular Hot body in the middle
        end
    end
end
Tinit=T0; % Store the initial conditions

timestep=0.75*((xstp^2+ystp^2)/2.0)/kappa/3.0; % Explicit timestep

% We have to create now a velocity field. Let’s try with a rotation
% around the center of the computational domain:
Vel=1.e-9; %m/s
% Also the advection scheme has a stability criteria.
% Beyond this point, it becomes unstable. Maxstep can be set to 0.05
% or 0.75 in order to see the difference.
maxStep=0.75;

if Vel*timestep>maxStep*xstp
    timestep=maxStep*xstp/Vel;
end

% Defining rotation velocity structure
for xn = 1:1:xnum
    for yn = 1:1:ynum
        VX(xn,yn)= Vel*(yn-(ynum-1)/2)/((ynum-1)/2);
        VY(xn,yn)= -Vel*(xn-(xnum-1)/2)/((xnum-1)/2);
    end
end