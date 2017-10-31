clear all
close all

syms c
%Moelans J/m³ 

% ABC_cu=[1e8 0 0];
% ABC_eta=[1e9 0 -1e6];
% ABC_sn=[1e9 0 1e7];
% 
% xhat_cu=0.076;
% xhat_eta=0.455;
% xhat_sn=0.978;

% %Min J/mol rumstemp
ABC_cu=[1.0133e5 -2.1146e4 -1.2842e4+1.9185e4]/16.29e-6;
ABC_eps=[0 0 0];
ABC_eta=[4e5 -6.9892e3 -1.9185e4+1.9185e4]/16.29e-6;
ABC_sn=[4.2059e6 7.1680e3 -1.5265e4+1.9185e4]/16.29e-6;

xhat_cu=0.10569;
xhat_eps = 0;
xhat_eta=0.41753;
xhat_sn=0.99941;

% 140 Li2009
% ABC_cu=[4.0359e5 -4.1474e4 -1.468e4];
% ABC_eta=[4e6 1.9832e3 -2.4536e4];
% ABC_sn=[7.3613e5 5.0136e3 -2.173e4];
% 
% xhat_cu=0.0109;
% xhat_eta=0.4355;
% xhat_sn=0.9953;

%Min J/m^3 220
% ABC_cu=[1.9926e10 -2.69e9 -1.1364e9];
% ABC_eps=[1.2277e11 2.2825e10 3.8822e8];
% ABC_eta=[2.4555e10 2.2279e7 -1.7646e9];
% ABC_sn=[2.3172e10 2.1576e8 -1.646e9];
% 
% xhat_cu=0.0171;
% xhat_eps=0.3281;
% xhat_eta=0.4359;
% xhat_sn=0.989;

G_cu=0.5*ABC_cu(1)*(c-xhat_cu).^2+ABC_cu(2)*(c-xhat_cu)+ABC_cu(3);
G_eps=0.5*ABC_eps(1)*(c-xhat_eps).^2+ABC_eps(2)*(c-xhat_eps)+ABC_eps(3);
G_eta=0.5*ABC_eta(1)*(c-xhat_eta).^2+ABC_eta(2)*(c-xhat_eta)+ABC_eta(3);
G_sn=0.5*ABC_sn(1)*(c-xhat_sn).^2+ABC_sn(2)*(c-xhat_sn)+ABC_sn(3);
dG_cu=ABC_cu(1)*(c-xhat_cu)+ABC_cu(2);
dG_eps=ABC_eps(1)*(c-xhat_eps)+ABC_eps(2);
dG_eta=ABC_eta(1)*(c-xhat_eta)+ABC_eta(2);
dG_sn=ABC_sn(1)*(c-xhat_sn)+ABC_sn(2);

%Find common tangents
f_cu_eta=@(x)commontangent(x,ABC_cu,ABC_eta,xhat_cu,xhat_eta);
f_eta_sn=@(x)commontangent(x,ABC_eta,ABC_sn,xhat_eta,xhat_sn);
f_cu_sn=@(x)commontangent(x,ABC_cu,ABC_sn,xhat_cu,xhat_sn);
%Cu-Cu6Sn5
x0=[0,0.5];
x_cu_eta=fsolve(f_cu_eta,x0)
y_cu_eta=subs(dG_cu,x_cu_eta(1))*c+subs(G_cu,x_cu_eta(1))-subs(dG_cu,x_cu_eta(1))*x_cu_eta(1); 

%Cu6Sn5-Sn
x0=[0.4,1];
x_eta_sn=fsolve(f_eta_sn,x0)
y_eta_sn=subs(dG_eta,x_eta_sn(1))*c+subs(G_eta,x_eta_sn(1))-subs(dG_eta,x_eta_sn(1))*x_eta_sn(1); 

%Cu-Sn
x0=[0,1];
x_cu_sn=fsolve(f_cu_sn,x0)
y_cu_sn=subs(dG_cu,x_cu_sn(1))*c+subs(G_cu,x_cu_sn(1))-subs(dG_cu,x_cu_sn(1))*x_cu_sn(1); 


hold on

fplot(y_cu_eta,[0,1], 'k')
fplot(y_eta_sn,[0,1],'k')
fplot(y_cu_sn,[0,1],'k')
fplot(G_cu,[0,1])
fplot(G_eps,[0,1])
fplot(G_eta,[0,1])
fplot(G_sn,[0,1])

title('Gibbs energy')


