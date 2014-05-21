function [T,VX,VY] = solving_advol( T_in )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute advention diffusion equation for the displacement, volume 
%%  deformation
%%  It is solved using upwind scheme
%%
%%  Inputs:  1. T0 -> the initial volume
%%  Outputs: 1. T1 -> the deformed volume   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_iter = 200; % Defining the number of maximum iterations
[r c s] = size(T_in);

% Defining grid in the horizontal direction
xmin = 0.0; xmax = c-1; xnum = c;
ymin = 0.0; ymax = r-1; ynum = r;
zmin = 0.0; zmax = s-1; znum = s;

xstp = (xmax - xmin) / (xnum - 1);
ystp = (ymax - ymin) / (ynum - 1);
zstp = (zmax - zmin) / (znum - 1);

Tinit = T_in'; % Store the initial conditions

timestep = 0.001;

%% Initialize the displacement
c = .5;

[VX, VY, VZ] = computing_displacement3(xnum, ynum, znum, xmin:xstp:xmax, ymin:ystp:ymax, zmin:zstp:zmax, xstp, ystp, zstp);


%%%%%%%%%% Solving using upwind scheme

T0 = Tinit; % Reinitializes the vector that contains Temperatures
T1 = zeros(size(T0));

% Solving temperature equation by explicit method + UPWIND scheme
for iter = 0:max_iter
    for xn = 3:xnum-2
        for yn = 3:ynum-2
            for zn = 3:znum-2

            T1(xn,yn,zn) = T0(xn,yn,zn);
            
            vx_pl = max(VX(xn,yn,zn),0);
            vx_mi = min(VX(xn,yn,zn),0);
            
            vy_pl = max(VY(xn,yn,zn),0);
            vy_mi = min(VY(xn,yn,zn),0);

            vz_pl = max(VZ(xn,yn,zn),0);
            vz_mi = min(VZ(xn,yn,zn),0);
            
            Tx_mi = (2 * T0(xn+1,yn,zn) + 3 * T0(xn,yn,zn)   - 6 * T0(xn-1,yn,zn) +     T0(xn-2,yn,zn)) / (6 * xstp);
            Tx_pl = (-   T0(xn+2,yn,zn) + 6 * T0(xn+1,yn,zn) - 3 * T0(xn,yn,zn)   - 2 * T0(xn-1,yn,zn)) / (6 * xstp);
            
            Ty_mi = (2 * T0(xn,yn+1,zn) + 3 * T0(xn,yn,zn)   - 6 * T0(xn,yn-1,zn) +     T0(xn,yn-2,zn)) / (6 * ystp);
            Ty_pl = (-   T0(xn,yn+2,zn) + 6 * T0(xn,yn+1,zn) - 3 * T0(xn,yn,zn)   - 2 * T0(xn,yn-1,zn)) / (6 * ystp);

            Tz_mi = (2 * T0(xn,yn,zn+1) + 3 * T0(xn,yn,zn)   - 6 * T0(xn,yn,zn-1) +     T0(xn,yn,zn-2)) / (6 * zstp);
            Tz_pl = (-   T0(xn,yn,zn+2) + 6 * T0(xn,yn,zn+1) - 3 * T0(xn,yn,zn)   - 2 * T0(xn,yn,zn-1)) / (6 * zstp);
            
            % And now we add the new UPWIND advection in time of the Temperature field:
            % Advective upwind scheme in X direction
            T1(xn,yn,zn) = T1(xn,yn,zn) - timestep*(vx_pl * Tx_mi + vx_mi * Tx_pl);
            T1(xn,yn,zn) = T1(xn,yn,zn) - timestep*(vy_pl * Ty_mi + vy_mi * Ty_pl);
            T1(xn,yn,zn) = T1(xn,yn,zn) - timestep*(vz_pl * Tz_mi + vz_mi * Tz_pl);

            end
        end
    end
    T0 = T1;% Reloading solutions to T0
end
T = T1';


