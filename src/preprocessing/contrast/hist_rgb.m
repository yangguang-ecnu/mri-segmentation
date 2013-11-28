function [h x] = hist_rgb(I)

[rows cols ch] = size(I);
h = zeros(256,ch);

[h(:,1),x] = histdouble(I(:,:,1));
h(:,2) = histdouble(I(:,:,2));
h(:,3) = histdouble(I(:,:,3));



% figure;
% bar(g);
% x = 0:255; 
% stem(g(:,1),'r','Marker','none');
% hold on
% stem(g(:,2),'Marker','none');
% stem(g(:,3),'g','Marker','none');
% hold off

%axis([1 3 0 255 0 max(max(g))])
% plot(x,h(:,1),'r');
% hold on
% plot(x,h(:,2),'g');
% plot(x,h(:,3),'b');
% hold off
% hleg = legend('Red','Blue','Green');
% set(hleg,'Location','NorthWest')
% set(hleg,'Interpreter','none')