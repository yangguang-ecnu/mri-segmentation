function [T,VX,VY] = solving_ad3( T_in )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute advention diffusion equation for the displacement, image 
%%  deformation
%%  It is solved using upwind scheme
%%
%%  Inputs:  1. T0 -> the initial image
%%  Outputs: 1. T1 -> the deformed image   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_iter = 200; % Defining the number of maximum iterations
[r c] = size(T_in);

% Defining grid in the horizontal direction
xmin = 0.0; xmax = c-1; xnum = c;
ymin = 0.0; ymax = r-1; ynum = r;

xstp = (xmax - xmin) / (xnum - 1);
ystp = (ymax - ymin) / (ynum - 1);

Tinit = T_in'; % Store the initial conditions

timestep = 0.001;

%% Initialize the displacement
c = .5;

[VX, VY] = computing_displacement(xnum, ynum, xmin:xstp:xmax, ymin:ystp:ymax, xstp, ystp);


%%%%%%%%%% Solving using upwind scheme

T0 = Tinit; % Reinitializes the vector that contains Temperatures
T1 = zeros(size(T0));
% Solving temperature equation by explicit method + UPWIND scheme
for iter = 0:max_iter
    for xn = 3:xnum-2
        for yn = 3:ynum-2
            
            T1(xn,yn) = T0(xn,yn);
            
            vx_pl = max(VX(xn,yn),0);
            vx_mi = min(VX(xn,yn),0);
            
            vy_pl = max(VY(xn,yn),0);
            vy_mi = min(VY(xn,yn),0);
            
            Tx_mi = (2 * T0(xn+1,yn) + 3 * T0(xn,yn)   - 6 * T0(xn-1,yn) + T0(xn-2,yn)) / (6 * xstp);
            Tx_pl = (-   T0(xn+2,yn) + 6 * T0(xn+1,yn) - 3 * T0(xn,yn)   - 2 * T0(xn-1,yn)) / (6 * xstp);
            
            Ty_mi = (2 * T0(xn,yn+1) + 3 * T0(xn,yn)   - 6 * T0(xn,yn-1) + T0(xn,yn-2)) / (6 * ystp);
            Ty_pl = (-   T0(xn,yn+2) + 6 * T0(xn,yn+1) - 3 * T0(xn,yn)   - 2 * T0(xn,yn-1)) / (6 * ystp);
            
            % And now we add the new UPWIND advection in time of the Temperature field:
            % Advective upwind scheme in X direction
            T1(xn,yn) = T1(xn,yn) - timestep*(vx_pl * Tx_mi + vx_mi * Tx_pl);
            T1(xn,yn) = T1(xn,yn) - timestep*(vy_pl * Ty_mi + vy_mi * Ty_pl);
            
        end
    end
    T0 = T1;% Reloading solutions to T0
end
T = T1';
% figure;
% subplot(131);imshow(Tinit,[]);
% subplot(132);imshow(T1,[]);
% subplot(133);imshow(Tinit - T1,[]);

% plotting_ad;