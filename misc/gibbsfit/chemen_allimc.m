close all
clear all

T=25+273.15; %Temperature
R=8.3144621; %Gas constant

Vm=16.29e-6; %m^3/mole
Vmcu=7.124e-6;%7.11e-6;%Vm;
Vmeps=8.6e-6;%Vm;
Vmeta=10.6e-6;%Vm;
Vmsn=Vm;


version=2; %1 for Li2009, 2 for Shim1996

if version==1
    %% Thermodynamic data from Li, Du et al. 2009
    GHSER_Cu=-7770.458+130.485235*T-24.112392*T*log(T)-2.65684e-3*T^2+0.129223e-6*T^3+52478*T^-1;
    if T<505
        GHSER_Sn=-5855.135+65.443315*T-15.961*T*log(T)-0.0188702*T^2+3.121167e-6*T^3-61960*T^-1;
    elseif T<800
        GHSER_Sn=2524.724+4.005269*T-8.2590486*T*log(T)-16.814429e-3*T^2+2.623131e-6*T^3-1081244*T^-1-1.2307e25*T^-9;
    end
    
    %Liquid
    Gliq_Cu=5194.277+120.973331*T-24.112392*T*log(T)-0.00265684*T^2+1.29223e-7*T^3+52478*T^-1-5.8489e-21*T^7;
    if T<505
        Gliq_Sn=1247.957+51.355548*T-15.961*T*log(T)-18.8702e-3*T^2+3.121167e-6*T^3-61960*T^-1+147.031e-20*T^7;
    elseif T<800
        Gliq_Sn=9496.31-9.809114*T-8.2590486*T*log(T)-16.814429e-3*T^2+2.623131e-6*T^3-1081244*T^-1;
    end
    L0liq=-8266.6-6.9973*T;
    L1liq=-21662+8.4655*T;
    L2liq=-26957.2+12.8887*T;
    L3liq=-13222.7+10.0420*T;
    
    %fcc
    Gfcc_Cu=GHSER_Cu;
    if T<505
        Gfcc_Sn=-345.135+56.983315*T-15.961*T*log(T)-0.0188702*T^2+3.121167e-6*T^3-61960*T^-1;
    elseif T<800
        Gfcc_Sn=8034.724-4.454731*T-8.2590486*T*log(T)-0.016814429*T^2+2.623131e-6*T^3-1081244*T^-1-1.2307e25*T^-9;
    end
    L0fcc=-11215.8+1.0978*T;
    L1fcc=-12884.1+5.4140*T;
    
    %bct
    Gbct_Cu=GHSER_Cu+5888;
    Gbct_Sn=GHSER_Sn;
    
elseif version==2
    %% Thermodynamic data from Shim1996
    GHSER_Cu=-7770.458+130.485235*T-24.112392*T*log(T)-2.65684e-3*T^2+0.129223e-6*T^3+52478*T^-1;
    if T<505
        GHSER_Sn=-5855.135+65.443315*T-15.961*T*log(T)-0.0188702*T^2+3.121167e-6*T^3-61960*T^-1;
    elseif T<800
        GHSER_Sn=2524.724+4.005269*T-8.2590486*T*log(T)-16.814429e-3*T^2+2.623131e-6*T^3-1081244*T^-1-1.2307e25*T^-9;
    end
    
    %Liquid
    Gliq_Cu=5194.277+120.973331*T-24.112392*T*log(T)-0.00265684*T^2+1.29223e-7*T^3+52478*T^-1-5.8489e-21*T^7;
    if T<505
        Gliq_Sn=1247.957+51.355548*T-15.961*T*log(T)-18.8702e-3*T^2+3.121167e-6*T^3-61960*T^-1+147.031e-20*T^7;
    elseif T<800
        Gliq_Sn=9496.31-9.809114*T-8.2590486*T*log(T)-16.814429e-3*T^2+2.623131e-6*T^3-1081244*T^-1;
    end
    L0liq=-9002.8-5.8381*T;
    L1liq=-20100.4+3.6366*T;
    L2liq=-10528.4;
    L3liq=0;
    
    %fcc
    Gfcc_Cu=GHSER_Cu;
    if T<505
        Gfcc_Sn=-1705.135+60.243315*T-15.961*T*log(T)-18.8702e-3*T^2+3.121167e-6*T^3-61960*T^-1;
    elseif T<800
        Gfcc_Sn=6674.724-1.194731*T-8.2590486*T*log(T)-16.814429e-3*T^2+2.623131e-6*T^3-1081244*T^-1-1.2307e25*T^-9;
    end
    L0fcc=-10672-1.4837*T;
    L1fcc=-15331.3+6.9539*T;
    
    %bct
    Gbct_Cu=0;%GHSER_Cu+5888;
    Gbct_Sn=GHSER_Sn;
end
%% Gibbs energy curves
syms x x1 x2

logterm=R*T*((1-x)*log(1-x)+x*log(x));

%Cu (fcc)
G_Cu=(1-x)*Gfcc_Cu+x*Gfcc_Sn+logterm+x*(1-x)*(L0fcc+L1fcc*(1-2*x));
% G_Cu=G_Cu/Vmcu;
%Sn (liq)
G_Snliq=(1-x)*Gliq_Cu+x*Gliq_Sn+logterm+x*(1-x)*(L0liq+L1liq*(1-2*x)+L2liq*(1-2*x)^2+L3liq*(1-2*x)^3);
% G_Snliq=G_Snliq/Vmsn;
%Sn (bct)
G_Snbct=(1-x)*Gbct_Cu+x*Gbct_Sn+logterm;
% G_Snbct=G_Snbct/Vmsn;
%As in Park Arroyave but with data from Du
% if version==1
%     G_eps=2e5*(x-0.249)^2+3*GHSER_Cu+GHSER_Sn-32690.5-0.7960*T;
%     G_eta=2e5*(x-0.435)^2+2*GHSER_Cu+GHSER_Sn-11566.8+8.8340*T;
% elseif version==2
    % %IMC from Park Arroyave
    %Cu3Sn
    G_eps=2e6*(x-0.249)^2+0.75*Gfcc_Cu+0.25*Gbct_Sn-8194.2-0.2043*T;
%     G_eps=G_eps/Vmeps;
    % %Cu6Sn5
    G_eta=2e6*(x-0.435)^2+0.545*Gfcc_Cu+0.455*Gbct_Sn-6869.5-0.1589*T;
%     G_eta=G_eta/Vmeta;
% end



% Derivatives
%Cu
mu_Cu=diff(G_Cu);
d2_Cu=diff(G_Cu,2);

%Sn (liq)
mu_Snliq=diff(G_Snliq);
d2_Snliq=diff(G_Snliq,2);

%Sn (bct)
mu_Snbct=diff(G_Snbct);
d2_Snbct=diff(G_Snbct,2);

%Cu3Sn
mu_eps=diff(G_eps);
d2_eps=diff(G_eps,2);

%Cu6Sn5
mu_eta=diff(G_eta);
d2_eta=diff(G_eta,2);

%% Equilibrium molar fractions
% %Cu-Sn liq
% init_guess=[0,1;0,1];
% e1=subs(G_Cu,x,x1)-subs(mu_Cu,x,x1)*x1==subs(G_Snliq,x,x2)-subs(mu_Snliq,x,x2)*x2; %Points on the same line
% e2=subs(mu_Cu,x,x1)==subs(mu_Snliq,x,x2); %with same derivative
% [x_cusnliq,x_snliqcu]=vpasolve([e1,e2],[x1,x2],init_guess);
% t1=subs(G_Cu,x,x_cusnliq)-subs(mu_Cu,x,x_cusnliq)*x_cusnliq+subs(mu_Cu,x,x_cusnliq)*x;

%Cu-Sn bct
init_guess=[0,1;0,1];
e1=subs(G_Cu,x,x1)-subs(mu_Cu,x,x1)*x1==subs(G_Snbct,x,x2)-subs(mu_Snbct,x,x2)*x2; %Points on the same line
e2=subs(mu_Cu,x,x1)==subs(mu_Snbct,x,x2); %with same derivative
[x_cusnbct,x_snbctcu]=vpasolve([e1,e2],[x1,x2],init_guess);
t2=subs(G_Cu,x,x_cusnbct)-subs(mu_Cu,x,x_cusnbct)*x_cusnbct+subs(mu_Cu,x,x_cusnbct)*x;

%Cu-Cu3Sn
init_guess=[0,0.5;0,0.5];
e1=subs(G_Cu,x,x1)-subs(mu_Cu,x,x1)*x1==subs(G_eps,x,x2)-subs(mu_eps,x,x2)*x2; %Points on the same line
e2=subs(mu_Cu,x,x1)==subs(mu_eps,x,x2); %with same derivative
[x_cueps,x_epscu]=vpasolve([e1,e2],[x1,x2],init_guess);
t3=subs(G_Cu,x,x_cueps)-subs(mu_Cu,x,x_cueps)*x_cueps+subs(mu_Cu,x,x_cueps)*x;

%Cu-Cu6Sn5
init_guess=[0,1;0,0.6];
e1=subs(G_Cu,x,x1)-subs(mu_Cu,x,x1)*x1==subs(G_eta,x,x2)-subs(mu_eta,x,x2)*x2; %Points on the same line
e2=subs(mu_Cu,x,x1)==subs(mu_eta,x,x2); %with same derivative
[x_cueta,x_etacu]=vpasolve([e1,e2],[x1,x2],init_guess);
t4=subs(G_Cu,x,x_cueta)-subs(mu_Cu,x,x_cueta)*x_cueta+subs(mu_Cu,x,x_cueta)*x;


%Cu3Sn-Cu6Sn5
init_guess=[0,1;0.3,1];
e1=subs(G_eps,x,x1)-subs(mu_eps,x,x1)*x1==subs(G_eta,x,x2)-subs(mu_eta,x,x2)*x2; %Points on the same line
e2=subs(mu_eps,x,x1)==subs(mu_eta,x,x2); %with same derivative
[x_epseta,x_etaeps]=vpasolve([e1,e2],[x1,x2],init_guess);
t5=subs(G_eps,x,x_epseta)-subs(mu_eps,x,x_epseta)*x_epseta+subs(mu_eps,x,x_epseta)*x;

% %Cu6Sn5-Sn liq
% init_guess=[0.5,1;0.3,1];
% e1=subs(G_Snliq,x,x1)-subs(mu_Snliq,x,x1)*x1==subs(G_eta,x,x2)-subs(mu_eta,x,x2)*x2; %Points on the same line
% e2=subs(mu_Snliq,x,x1)==subs(mu_eta,x,x2); %with same derivative
% [x_snliqeta,x_etasnliq]=vpasolve([e1,e2],[x1,x2],init_guess);
% t6=subs(G_Snliq,x,x_snliqeta)-subs(mu_Snliq,x,x_snliqeta)*x_snliqeta+subs(mu_Snliq,x,x_snliqeta)*x;

%Cu6Sn5-Sn bct
init_guess=[0.5,1;0.3,1];
e1=subs(G_Snbct,x,x1)-subs(mu_Snbct,x,x1)*x1==subs(G_eta,x,x2)-subs(mu_eta,x,x2)*x2; %Points on the same line
e2=subs(mu_Snbct,x,x1)==subs(mu_eta,x,x2); %with same derivative
[x_snbcteta,x_etasnbct]=vpasolve([e1,e2],[x1,x2],init_guess);
t7=subs(G_Snbct,x,x_snbcteta)-subs(mu_Snbct,x,x_snbcteta)*x_snbcteta+subs(mu_Snbct,x,x_snbcteta)*x;

%% Fitting
%Cu
% xhat_Cusnliq=(x_cusnliq+x_cueps)/2;
% A_Cusnliq=subs(d2_Cu,x,xhat_Cusnliq);
% B_Cusnliq=subs(mu_Cu,x,xhat_Cusnliq);
% C_Cusnliq=subs(G_Cu,x,xhat_Cusnliq);
%
% Gf_Cusnliq=A_Cusnliq*(x-xhat_Cusnliq)^2+B_Cusnliq*(x-xhat_Cusnliq)+C_Cusnliq;

xhat_Cusnbct=x_cueps;%(x_cusnbct+x_cueps)/2;
A_Cusnbct=subs(d2_Cu,x,xhat_Cusnbct);
B_Cusnbct=subs(mu_Cu,x,xhat_Cusnbct);
C_Cusnbct=subs(G_Cu,x,xhat_Cusnbct);

Gf_Cusnbct=0.5*A_Cusnbct*(x-xhat_Cusnbct)^2+B_Cusnbct*(x-xhat_Cusnbct)+C_Cusnbct;

% % Cu3Sn
xhat_eps=0.5*(x_epscu+x_epseta);
A_eps=subs(d2_eps,x,xhat_eps);
B_eps=subs(mu_eps,x,xhat_eps);
C_eps=subs(G_eps,x,xhat_eps);

Gf_eps=0.5*A_eps*(x-xhat_eps)^2+B_eps*(x-xhat_eps)+C_eps;

% Cu6Sn5
% xhat_etasnliq=(x_etasnliq+x_etaeps)/2;
% A_etasnliq=subs(d2_eta,x,xhat_etasnliq);
% B_etasnliq=subs(mu_eta,x,xhat_etasnliq);
% C_etasnliq=subs(G_eta,x,xhat_etasnliq);
%
% Gf_etasnliq=0.5*A_etasnliq*(x-xhat_etasnliq)^2+B_etasnliq*(x-xhat_etasnliq)+C_etasnliq;

xhat_etasnbct=(x_etasnbct+x_etaeps)/2;
A_etasnbct=subs(d2_eta,x,xhat_etasnbct);
B_etasnbct=subs(mu_eta,x,xhat_etasnbct);
C_etasnbct=subs(G_eta,x,xhat_etasnbct);

Gf_etasnbct=0.5*A_etasnbct*(x-xhat_etasnbct)^2+B_etasnbct*(x-xhat_etasnbct)+C_etasnbct;

% Sn
% xhat_snliq=(x_snliqeta+x_snliqcu)/2;
% A_Snliq=subs(d2_Snliq,x,xhat_snliq);
% B_Snliq=subs(mu_Snliq,x,xhat_snliq);
% C_Snliq=subs(G_Snliq,x,xhat_snliq);
%
% Gf_Snliq=0.5*A_Snliq*(x-xhat_snliq)^2+B_Snliq*(x-xhat_snliq)+C_Snliq;


xhat_snbct=x_snbcteta;%(x_snbcteta+x_snbctcu)/2;
A_Snbct=subs(d2_Snbct,x,xhat_snbct);
B_Snbct=subs(mu_Snbct,x,xhat_snbct);
C_Snbct=subs(G_Snbct,x,xhat_snbct);

Gf_Snbct=0.5*A_Snbct*(x-xhat_snbct)^2+B_Snbct*(x-xhat_snbct)+C_Snbct;

%% Plotting
figure()
hold on
p1=ezplot(G_Cu,[0,1]);
p4=ezplot(G_Snbct,[0,1]);
p5=ezplot(Gf_Cusnbct,[0,1]);
p6=ezplot(Gf_eps,[0,1]);
p7=ezplot(Gf_etasnbct,[0,1]);
p8=ezplot(Gf_Snbct,[0,1]);
title(' ')
ylabel('Gibbs energy [J/mol]')
xlabel('Mole fraction of Sn')
axis([0,1,-31e3,-12e3]);
legend('Cu','Sn bct','Cu fit','Cu3Sn fit','Cu6Sn5 fit','Sn bct fit','Location','se')
hold off



figure()
hold on
p1=ezplot(Gf_Cusnbct/Vmcu,[0,1]);
p2=ezplot(Gf_eps/Vmeps,[0,1]);
p3=ezplot(Gf_etasnbct/Vmeta,[0,1]);
p4=ezplot(Gf_Snbct/Vmsn,[0,1]);
p5=ezplot(G_Cu/Vmcu,[0,1]);
p8=ezplot(G_Snbct/Vmsn,[0,1]);
title('Sn bct')
ylabel('Gibbs energy [J/m^3]')
xlabel('Mole fraction of Sn')
axis([0,1,-3e9,0]);
legend('Cu-fit','Cu3Sn-fit','Cu6Sn5-fit','Sn bct-fit','Cu','Sn bct','Location','se')
hold off

%% Display
% disp('Cu snliq')
% A_Cusnliq=eval(A_Cusnliq)
% B_Cusnliq=eval(B_Cusnliq)
% C_Cusnliq=eval(C_Cusnliq)
% disp('--------')
disp('Cu snbct')
A_Cusnbct=eval(A_Cusnbct)
B_Cusnbct=eval(B_Cusnbct)
C_Cusnbct=eval(C_Cusnbct)
disp('--------')
disp('Cu3Sn snbct')
A_eps=eval(A_eps)
B_eps=eval(B_eps)
C_eps=eval(C_eps)
disp('--------')
% disp('Cu6Sn5 sn liq')
% A_etasnliq=eval(A_etasnliq)
% B_etasnliq=eval(B_etasnliq)
% C_etasnliq=eval(C_etasnliq)
% disp('---------')
disp('Cu6Sn5 sn bct')
A_etasnbct=eval(A_etasnbct)
B_etasnbct=eval(B_etasnbct)
C_etasnbct=eval(C_etasnbct)
% disp('---------')
% disp('Sn liq')
% A_Snliq=eval(A_Snliq)
% B_Snliq=eval(B_Snliq)
% C_Snliq=eval(C_Snliq)
disp('---------')
disp('Sn bct')
A_Snbct=eval(A_Snbct)
B_Snbct=eval(B_Snbct)
C_Snbct=eval(C_Snbct)
disp('-----------')
disp('Interface concentrations')
% xhat_Cusnliq=eval(xhat_Cusnliq)
xhat_Cusnbct=eval(xhat_Cusnbct)
xhat_eps=eval(xhat_eps)
% xhat_etasnliq=eval(xhat_etasnliq)
xhat_etasnbct=eval(xhat_etasnbct)
% xhat_Snliq=eval(xhat_snliq)
xhat_Snbct=eval(xhat_snbct)

xcu_eps=eval(x_cueps)
xcu_eta=eval(x_cueta)
xcu_sn=eval(x_cusnbct)
xeps_cu=eval(x_epscu)
xeps_eta=eval(x_epseta)
xeta_cu=eval(x_etacu)
xeta_eps=eval(x_etaeps)
xeta_sn=eval(x_etasnbct)
xsn_cu=eval(x_snbctcu)
xsn_eta=eval(x_snbcteta)