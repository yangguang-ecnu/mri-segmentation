% % Define and fill a rectangular area in the plane
% xlimit = [3 13];
% ylimit = [2  8];
% zlimit = [2  9];
% xbox = xlimit([1 1 2 2 1]);
% ybox = ylimit([1 2 2 1 1]);
% zbox = zlimit([2 2 1 1 1]);
% %mapshow(xbox,ybox,zbox,'DisplayType','polygon','LineStyle','none')
% 
% % Define and display a two-part polyline
% x = [0 6  4  8 8 10 14 10 14 NaN 4 4 6 9 15];
% y = [4 6 10 11 7  6 10 10  6 NaN 0 3 4 3  6];
% z = [4 6 10 11 7  6 10 10  6 NaN 0 3 4 3  6];
% %mapshow(x,y,'Marker','+')
% 
% % Intersect the polyline with the rectangle
% [xi, yi, zi] = polyxpoly(x, y, z, xbox, ybox, zbox);
% %mapshow(xi,yi,'DisplayType','point','Marker','o')

function PolygonClip_example

    P1.x=[-1 1 1 -1]; P1.y=[-1 -1 1 1]; P1.z = [1 1 1 1]; P1.hole=0;
    P1(2).x=[-1 1 1 -1]*.5; P1(2).y=[-1 -1 1 1]*.5; P1(2).z = [1 1 1 1];P1(2).hole=1;

    P2.x=[-2 0.8 0.4 -2]; P2.y=[-.5 -.2 0 .5]; P2.z = [0 0 0 0];P2.hole=0;
    P2(2).x=[2 0.8 0.6 1.5]; P2(2).y=[-1 0 0.3 1]; P2(2).z = [0 0 0 0];P2(2).hole=0;

    for type = 0:3
        subplot(2,2,type+1); box on
        switch type
            case 0; title('A-B')
            case 1; title('A.and.B (standard)')
            case 2; title('xor(A,B)')
            case 3; title('union(A,B)')
        end

        P3=PolygonClip(P1,P2,type);

        for i=1:3
            eval(['p=P' num2str(i) ';'])
            for np=1:length(p)
                obj=patch(p(np).x,p(np).y,i);
                if p(np).hole==1; set(obj,'facecolor','w'); end
            end
        end

    end