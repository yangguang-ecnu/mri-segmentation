function T1 = solving_ad( T0 )
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

max_iter = 100; % Defining the number of maximum iterations
[r c] = size(T0);

% Defining grid in the horizontal direction
xmin = 0.0; xmax = c-1; xnum = c;
ymin = 0.0; ymax = r-1; ynum = r;

xstp = (xmax - xmin) / (xnum - 1);
ystp = (ymax - ymin) / (ynum - 1);

Tinit = T0; % Store the initial conditions

timestep = 1e-2;

%% Initialize the displacement
c = .5;

[VX, VY] = computing_displacement(xnum, ynum, xmin:xstp:xmax, ymin:ystp:ymax, 1, 1, c);


%%%%%%%%%% Solving using upwind scheme

T0 = Tinit; % Reinitializes the vector that contains Temperatures

% Solving temperature equation by explicit method + UPWIND scheme
for iter = 0:max_iter
    for xn = 1:xnum
        for yn = 1:ynum
            % Boundary conditions:
            if (xn==1 || xn==xnum || yn==1 || yn==ynum)
                T1(xn,yn) = 0; % T=const
            else

                T1(xn,yn) = T0(xn,yn);
                
                % And now we add the new UPWIND advection in time of the Temperature field:
                % Advective upwind scheme in X direction
                if (VX(xn,yn)>0)
                    T1(xn,yn) = T1(xn,yn) - VX(xn,yn)*timestep*(T0(xn,yn)   - T0(xn-1,yn))/xstp;
                end
                if (VX(xn,yn)<0)
                    T1(xn,yn) = T1(xn,yn) - VX(xn,yn)*timestep*(T0(xn+1,yn) - T0(xn,yn))/xstp;
                end 
                % Advective upwind scheme in X direction
                if (VY(xn,yn)>0)
                    T1(xn,yn) = T1(xn,yn) - VY(xn,yn)*timestep*(T0(xn,yn)   - T0(xn,yn-1))/ystp;
                end
                if (VY(xn,yn)<0)
                    T1(xn,yn) = T1(xn,yn) - VY(xn,yn)*timestep*(T0(xn,yn+1) - T0(xn,yn))/ystp;
                end
            end
        end
    end
    T0 = T1;% Reloading solutions to T0
end

% figure;
% subplot(131);imshow(Tinit,[]);
% subplot(132);imshow(T1,[]);
% subplot(133);imshow(Tinit - T1,[]);
% 
% plotting_ad;