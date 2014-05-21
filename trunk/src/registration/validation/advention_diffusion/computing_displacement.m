function [u, v] = computing_displacement(nx, ny, x, y, dx, dy)

%Specifying Parameters

u = zeros(nx,ny);                  %Preallocating u
v = zeros(nx,ny);                  %Preallocating v

center = 0;%nx / 4;

% Initial Conditions
for i=1:nx
    for j=1:ny
        if ((center <= y(j))&&(y(j) <= ny/2)&&(center <= x(i))&&(x(i) <= nx/2))
            u(i,j)=27;
            v(i,j)=23;
        elseif ((ny >= y(j))&&(y(j) >= ny/2)&&(nx >= x(i))&&(x(i) >= nx/2))
            u(i,j)=17;
            v(i,j)=13;            
        else
            u(i,j)=10;
            v(i,j)=10;
        end
    end
end

% Boundary conditions
u(1,:)=0;
u(nx,:)=0;
u(:,1)=0;
u(:,ny)=0;
v(1,:)=0;
v(nx,:)=0; 
v(:,1)=0;
v(:,ny)=0;

i=2:nx-1;
j=2:ny-1;

dt = 0.01; 
nt = 100;
%%
%Explicit method with F.D in time and B.D in space
% figure;
for it=0:nt
    un=u;
    vn=v;
%     h=quiver(x(1:20:end,1:20:end),y(1:20:end,1:20:end),u(1:20:end,1:20:end)',v(1:20:end,1:20:end)','Color','black');       %plotting the velocity field
%     axis([0 nx 0 ny])
%     title({'2-D Convection';'Transport property vector field {\bfu}=(u_x,u_y)';['time(\itt) = ',num2str(dt*it)]})
%     xlabel('Spatial co-ordinate (x) \rightarrow')
%     ylabel('Spatial co-ordinate (y) \rightarrow')
%     drawnow; 
%     refreshdata(h)
    u(i,j) = un(i,j)-(dt*un(i,j).*(un(i,j) - un(i-1,j))/dx)-(dt*vn(i,j).*(un(i,j) - un(i,j-1))/dy);
    v(i,j) = vn(i,j)-(dt*un(i,j).*(vn(i,j) - vn(i-1,j))/dx)-(dt*vn(i,j).*(vn(i,j) - vn(i,j-1))/dy);
    %Boundary Conditions
    u(1,:)=0;
    u(nx,:)=0;
    u(:,1)=0;
    u(:,ny)=0;
    v(1,:)=0;
    v(nx,:)=0;
    v(:,1)=0;
    v(:,ny)=0;
end