function [u, v, w] = computing_displacement3(nx, ny, nz, x, y, z, dx, dy, dz)

%Specifying Parameters

u = zeros(nx,ny,nz);                  %Preallocating u
v = zeros(nx,ny,nz);                  %Preallocating v
w = zeros(nx,ny,nz);                  %Preallocating w

center = 0;%nx / 4;

% Initial Conditions
for i=1:nx
    for j=1:ny
        for k=1:nz
            if ((center <= y(j))&&(y(j) <= ny/2)&&(center <= x(i))&&(x(i) <= nx/2)&&(center <= z(k))&&(z(k) <= nz/2))
                u(i,j,k)=27;
                v(i,j,k)=23;
                v(i,j,k)=23;
            elseif ((ny >= y(j))&&(y(j) >= ny/2)&&(nx >= x(i))&&(x(i) >= nx/2)&&(nz >= z(k))&&(z(k) >= nz/2))
                u(i,j,k)=17;
                v(i,j,k)=13;
                v(i,j,k)=23;           
            else
                u(i,j,k)=10;
                v(i,j,k)=10;
                v(i,j,k)=23;
            end
        end
    end
end

% Boundary conditions
%% u
u(1, :,:)=0;
u(nx,:,:)=0;
u(:, 1,:)=0;
u(:,ny,:)=0;
u(:,:, 1)=0;
u(:,:,nz)=0;
%% v
v(1, :,:)=0;
v(nx,:,:)=0; 
v(:, 1,:)=0;
v(:,ny,:)=0;
v(:,:, 1)=0;
v(:,:,nz)=0;
%% w
w(1, :,:)=0;
w(nx,:,:)=0; 
w(:, 1,:)=0;
w(:,ny,:)=0;
w(:,:, 1)=0;
w(:,:,nz)=0;

i=2:nx-1;
j=2:ny-1;
k=2:nz-1;

dt = 0.01; 
nt = 100;

%Explicit method with F.D in time and B.D in space
for it=0:nt
    un = u;
    vn = v;
    wn = w;

    u(i,j,k) = un(i,j,k)-(dt*un(i,j,k).*(un(i,j,k) - un(i-1,j,k))/dx)-(dt*vn(i,j,k).*(un(i,j,k) - un(i,j-1,k))/dy)-(dt*wn(i,j,k).*(un(i,j,k) - un(i,j,k-1))/dz);
    v(i,j,k) = vn(i,j,k)-(dt*un(i,j,k).*(vn(i,j,k) - vn(i-1,j,k))/dx)-(dt*vn(i,j,k).*(vn(i,j,k) - vn(i,j-1,k))/dy)-(dt*wn(i,j,k).*(vn(i,j,k) - vn(i,j,k-1))/dz);
    w(i,j,k) = wn(i,j,k)-(dt*un(i,j,k).*(wn(i,j,k) - wn(i-1,j,k))/dx)-(dt*vn(i,j,k).*(wn(i,j,k) - wn(i,j-1,k))/dy)-(dt*wn(i,j,k).*(wn(i,j,k) - wn(i,j,k-1))/dz);

    % Boundary conditions
    %% u
    u(1, :,:)=0;u(nx,:,:)=0;u(:, 1,:)=0;u(:,ny,:)=0;u(:,:, 1)=0;u(:,:,nz)=0;
    %% v
    v(1, :,:)=0;v(nx,:,:)=0; v(:, 1,:)=0;v(:,ny,:)=0;v(:,:, 1)=0;v(:,:,nz)=0;
    %% w
    w(1, :,:)=0;w(nx,:,:)=0; w(:, 1,:)=0;w(:,ny,:)=0;w(:,:, 1)=0;w(:,:,nz)=0;
end


