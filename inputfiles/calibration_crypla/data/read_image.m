%--------------------------------------------------------------------------
%
% READ ME
% 
% 1. Pick three points on screen, the first one being the origin, x-axis,
% y-axis
% 2. Pick any number of points in, e.g., diagram
% 3. Press enter when ready
%--------------------------------------------------------------------------
clc
clear all
%strainscale=25.536945812807858/0.05;
%--------------------------------------------------------------------------
% read image
%--------------------------------------------------------------------------
[X,map]=imread('kiener2006fig4b.png','png');
%--------------------------------------------------------------------------
figure(1)
clf
% image(X)
imshow(X)
%--------------------------------------------------------------------------
% reaf image data
%--------------------------------------------------------------------------
[x_1,x_2] = ginput
%%
a(:,1) = x_1;
a(:,2) = x_2;

a(:,1) = a(:,1)-a(1,1);
a(:,2) = a(:,2)-a(1,2);
%--------------------------------------------------------------------------
% scale data points w.r.t. image axes
%--------------------------------------------------------------------------
picked_x_1 = input('Enter picked x_1 axis value: ');%-2;
picked_x_2 = input('Enter picked y_2 axis value: ');
scale_x_1 = a(2,1)/picked_x_1;
scale_x_2 = a(3,2)/picked_x_2;
a(:,1) = a(:,1)/scale_x_1;
a(:,2) = a(:,2)/scale_x_2;


%--------------------------------------------------------------------------
% plot result
%--------------------------------------------------------------------------
xx = a(4:end,1);
yy = a(4:end,2);


% %move data to origin
% xx=xx-xx(1);
% %remove "extra" point
% xx(1)=[];
% yy(1)=[];
figure()
plot(xx,yy)

% pp = interp1(a(4:end,1),a(4:end,2),'cubic','pp');
% yy = ppval(pp,xx);
% yy(end)=yy(1);
% yy(end-1)=yy(2);
% yy(end-2)=yy(3);
% plot(a(4:end,1),a(4:end,2),'o',xx,yy);
% save('pointssn.mat')
% save('pointscu.mat')