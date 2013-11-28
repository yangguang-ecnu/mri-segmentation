function plot_axes_p(axes_h,x,y,n)



if n==1
    axes(axes_h)
    plot(x,y);
    axis([0 max(x) 0 max(max(y))])
elseif n==3
    axes(axes_h)
    plot(x,y(:,1),'r');
    hold on
    plot(x,y(:,2),'g');
    plot(x,y(:,3),'b');
    hold off
    axis([0 max(x) 0 max(max(y))])
    hleg = legend('Red','Blue','Green');
    set(hleg,'Location','NorthWest')
    set(hleg,'Interpreter','none')
end