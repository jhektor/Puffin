%clear all
close all
%clc

% set options for optimization
option = optimset('Display', 'iter', 'FunValCheck', 'on',...
    'MaxIter',600, 'MaxFunEvals', 3000, 'TolFun', 1e-12, 'TolX', 1e-12,...
    'PlotFcns', {@optimplotfval, @optimplotx});

%load experimental values
lcase=2
% if lcase==1
%   S=load('../data/darbandi2012.mat','xx','yy');
%   expstrain=S.xx;
%   expstress=S.yy;
% elseif lcase==2 %philippi 110
%   S=load('../data/philippi2016extra110.mat','xx','yy');
%   expstrain=S.xx;
%   expstress=S.yy;
% elseif lcase==3 %philippi 001
%   S=load('../data/philippi2016extra001.mat','xx','yy');
%   expstrain=S.xx;
%   expstress=S.yy;
% elseif lcase==4 %Ekinci 110 fig1
%   S=load('../data/ekinci2003fig1-110.mat','xx','yy');
%   expstrain=S.xx;
%   expstress=1e6*S.yy;
% end
% starting values
if lcase==5 || lcase==6 %Cu
    E=110e9%93e9;%120e9;
    v=0.34;%0.35;
    
    G=E/(2*(1+v));   
    K=E/(3*(1-2*v));   
    
    G0=40e6;%64e6;%95e6;%27e6;
    gamma0=1e-3;%2e-5;%1.5e-5;%1e-4;
    mint=15;%16;%26
    hab=1.4;
    QQ=490e6;%261e6;%125e6;%100e6
    BB=15;%24;%14%8;
    gg=0.007;
else
    E=50e9;%50e9;
    v=0.36;
    
    G=E/(2*(1+v));   %18e9%6e9;%20e9;
    K=E/(3*(1-2*v));   %58e9;%25e9;%52e9;
    G0=27e6;%45e6;%50e6;%7e6;%5.6e6;%8.8e6;%9.3e6;%5e6;%5e8;
    gamma0=1e-3;%0.002;%0.0017;%1e-3;%1e-2;
    mint=8;%10;%8;%20;%5.;
    hab=1.4;
    QQ=26e6;%40e6;%2e6;%16e6;%20e6%2e6;%2.8e6;%1e4;
    BB=8;%20;%0.8;%0.7;%1;%0.01;
    gg=0.007;%0.03;%0.02;%0.01
    lb=[49e9,0.35,0,9.9e-4,8,1.39,0,7,0.006];
    ub=[51e9,0.37,Inf,1.1e-3,Inf,1.41,Inf,9,0.008];
    
    lb=zeros(9,1);
    ub=Inf*ones(9,1);
%     lb(1)=49e9; ub(1)=51e9; %E
    lb(2)=0.34; ub(2)=0.38; %v
%     lb(4)=9.9e-4; ub(4)=1.1e-3; %gamma0
    lb(6)=1.39; ub(6)=1.41; %hab
% %     lb(8)=7;    ub(8)=9; %B
%     lb(9)=0.006;ub(9)=0.008;
end
x0=[E,v,G0,gamma0,mint,hab,QQ,BB,gg];



%err=minimize_error(x0,expstrain,expstress)

[x,fval,~,options] = fminsearch(@(x) minimize_error(x,1,lb,ub), x0, option);
%[x,fval,exitflag,info] = fmincon(@(x) minimize_error(x,1), x0,[],[],[],[],lb,ub,[],option);
%x=x0
E=x(1)
v=x(2)
G=E/(2*(1+v));   %18e9%6e9;%20e9;
K=E/(3*(1-2*v));
% E=9*G*K/(3*K+G)
% v=(3*K-2*G)/(6*K+2*G)

G0=x(3)
gamma0=x(4)
mint=x(5)
hab=x(6)
QQ=x(7)
BB=x(8)
gg=x(9)
