function [Vertices Lines] = contour(position)
verbose = 1;
nBetween = 3;

p.x = position(:,2);
p.y = position(:,1);
p.n = size(position,1);

p.x(end +1) = p.x(1);
p.y(end +1) = p.y(1);

% Interpolate to get more points
r=10;
pointsx=interp(p.x,r); pointsx=pointsx(1:end-r+1);
pointsy=interp(p.y,r); pointsy=pointsy(1:end-r+1);


%imshow(handles.views.axial(:,:,handles.current_axial),[]);hold on
plot(pointsy,pointsx,'b','LineWidth',2); 

totalx=[]; totaly=[];
pointst = 1:length(pointsx);
i=find(pointst);
% Loop to make points evenly spaced on line pieces between landmark points
for j=1:length(i)-1
    % One line piece
    linex=pointsx(i(j):i(j+1));
    liney=pointsy(i(j):i(j+1));
    % Lenght on line through the points
    dx=[0 linex(2:end)-linex(1:end-1)];
    dy=[0 liney(2:end)-liney(1:end-1)];
    dist=cumsum(sqrt(dx.^2+dy.^2));
    % Interpolate new points evenly spaced on the line piece
    dist2=linspace(0,max(dist),nBetween);
    linex=interp1(dist,linex,dist2);
    liney=interp1(dist,liney,dist2);
    % Display the line piece
    %if(verbose),
        plot(liney,linex,'g*','MarkerSize',4);
        plot(liney(1),linex(1),'r*','MarkerSize',4);
        plot(liney(end),linex(end),'r*','MarkerSize',4); 
    %end
    % Remove Point because it is also in the next line piece
    if(j<i-1), linex(end)=[]; liney(end)=[]; end
    % Add the evenly spaced line piece to the total line
    totalx=[totalx linex];
    totaly=[totaly liney];
end
%Vertices = [totalx(:) totaly(:)];
Vertices = [p.y(:) p.x(:)];
Lines = [(1:size(Vertices,1))' ([2:size(Vertices,1) 1])'];