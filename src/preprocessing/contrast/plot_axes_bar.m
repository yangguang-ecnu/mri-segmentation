function plot_axes_bar(axes_h,x,y,n)



if n==1
    axes(axes_h)
    bar(x(1:3:end),y(1:3:end));
    axis([0 max(x) 0 max(max(y))])
elseif n==3
    axes(axes_h)
    bar(x(1:3:end),y(1:3:end,1),'r');
    hold on
    bar(x(1:3:end),y(1:3:end,2),'g');
    bar(x(1:3:end),y(1:3:end,3),'b');
    hold off
    axis([0 max(x) 0 max(max(y))])
    hleg = legend('Red','Blue','Green');
    set(hleg,'Location','NorthWest')
    set(hleg,'Interpreter','none')
end