function plot_axes_st(axes_h,x,y,n)



if n==1
    axes(axes_h)
    stem(x(1:3:end),y(1:3:end),'Marker','none');
    axis([0 max(x) 0 max(max(y))])
elseif n==3
    axes(axes_h)
    stem(x(1:3:end),y(1:3:end,1),'r','Marker','none');
    hold on
    stem(x(1:3:end),y(1:3:end,2),'g','Marker','none');
    stem(x(1:3:end),y(1:3:end,3),'b','Marker','none');
    hold off
    axis([0 max(x) 0 max(max(y))])
    hleg = legend('Red','Blue','Green');
    set(hleg,'Location','NorthWest')
    set(hleg,'Interpreter','none')
end