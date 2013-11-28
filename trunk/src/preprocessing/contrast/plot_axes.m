function plot_axes(axes_h,x,y,n,color)

disp(nargchk(1, 5, nargin));

if nargin == 4
    color = 'b';
end

if n==1
    axes(axes_h)
    stem(x,y,color,'Marker','none');
    axis([0 max(x) 0 max(max(y))])
elseif n==3
    axes(axes_h)
    stem(x,y(:,1),'r','Marker','none');
    hold on
    stem(x,y(:,2),'Marker','none');
    stem(x,y(:,3),'g','Marker','none');
    hold off
    axis([0 max(x) 0 max(max(y))])
    hleg = legend('Red','Blue','Green');
    set(hleg,'Location','NorthWest')
    set(hleg,'Interpreter','none')
end
