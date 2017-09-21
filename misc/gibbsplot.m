clear all
close all

syms c
ABC_cu=[1e8 0 0];
ABC_imc=[1e9 0 -1e6];
ABC_sn=[1e9 0 1e7];

xhat_cu=0.076;
xhat_imc=0.455;
xhat_sn=0.978;

G_cu=0.5*ABC_cu(1)*(c-xhat_cu).^2+ABC_cu(2)*(c-xhat_cu)+ABC_cu(3);
G_imc=0.5*ABC_imc(1)*(c-xhat_imc).^2+ABC_imc(2)*(c-xhat_imc)+ABC_imc(3);
G_sn=0.5*ABC_sn(1)*(c-xhat_sn).^2+ABC_sn(2)*(c-xhat_sn)+ABC_sn(3);
dG_cu=ABC_cu(1)*(c-xhat_cu)+ABC_cu(2);
dG_imc=ABC_imc(1)*(c-xhat_imc)+ABC_imc(2);
dG_sn=ABC_sn(1)*(c-xhat_sn)+ABC_sn(2);

%Find common tangents
f_cu_imc=@(x)commontangent(x,ABC_cu,ABC_imc,xhat_cu,xhat_imc);
f_imc_sn=@(x)commontangent(x,ABC_imc,ABC_sn,xhat_imc,xhat_sn);
f_cu_sn=@(x)commontangent(x,ABC_cu,ABC_sn,xhat_cu,xhat_sn);
%Cu-Cu6Sn5
x0=[0,0.5];
x_cu_imc=fsolve(f_cu_imc,x0)
y_cu_imc=subs(dG_cu,x_cu_imc(1))*c+subs(G_cu,x_cu_imc(1))-subs(dG_cu,x_cu_imc(1))*x_cu_imc(1); 

%Cu6Sn5-Sn
x0=[0.4,1];
x_imc_sn=fsolve(f_imc_sn,x0)
y_imc_sn=subs(dG_imc,x_imc_sn(1))*c+subs(G_imc,x_imc_sn(1))-subs(dG_imc,x_imc_sn(1))*x_imc_sn(1); 

%Cu-Sn
x0=[0,1];
x_cu_sn=fsolve(f_cu_sn,x0)
y_cu_sn=subs(dG_cu,x_cu_sn(1))*c+subs(G_cu,x_cu_sn(1))-subs(dG_cu,x_cu_sn(1))*x_cu_sn(1); 


hold on

fplot(y_cu_imc,[0,1], 'k')
fplot(y_imc_sn,[0,1],'k')
fplot(y_cu_sn,[0,1],'k')
fplot(G_cu,[0,1])
fplot(G_imc,[0,1])
fplot(G_sn,[0,1])

title('Gibbs energy')

