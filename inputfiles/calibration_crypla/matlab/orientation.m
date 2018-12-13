%Convert miller idices to Euler angles
% syms phi1 theta phi2 %alpha beta gamma
% R=[cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(theta), -cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(theta), sin(phi1)*sin(theta);
%    sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(theta), -sin(phi1)*sin(phi2)+cos(phi1)*cos(theta)*cos(phi2), -cos(phi1)*sin(theta);
%    sin(phi2)*sin(theta), cos(phi2)*sin(theta),cos(theta)];
% 
% hkl=[1,0,0]';
% eq1=R(1,:)*hkl==0;
% eq2=R(2,:)*hkl==0;
% eq3=R(3,:)*hkl==1;
% 
% 
% 
% Y=solve([eq1,eq2,eq3],[phi1,theta,phi2])
% 
% %Y=solve([R*hkl==[0;0;1]],[phi1,theta,phi2])

hkl=[1,1,0]'
hkl=hkl/norm(hkl);
x0=[90,45,0]*pi/180; %initial guess phi1,theta,phi2
f=@(x)[[cos(x(1))*cos(x(3))-sin(x(1))*sin(x(3))*cos(x(2)), -cos(x(1))*sin(x(3))-sin(x(1))*cos(x(3))*cos(x(2)), sin(x(1))*sin(x(2))]*hkl; ...
        [sin(x(1))*cos(x(3))+cos(x(1))*sin(x(3))*cos(x(2)), -sin(x(1))*sin(x(3))+cos(x(1))*cos(x(2))*cos(x(3)), -cos(x(1))*sin(x(2))]*hkl; ...
        [sin(x(3))*sin(x(2)), cos(x(3))*sin(x(2)),cos(x(2))]*hkl-1];

 opt=optimoptions('fsolve','maxfunevals',1000,'tolfun',1e-8)
    
[x,fval]=fsolve(f,x0,opt)

xdeg=180*x/pi