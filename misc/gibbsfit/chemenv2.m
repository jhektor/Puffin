%% This is the one used in the paper
%% Fit of Gibbs energy, as in Moelans 2015
close all
clear all

T=220+273.15;%453;%25+273.15;
R=8.3144621;

Vm=16.29e-6; %m^3/mole
Vmcu=Vm;%7.124e-6;%7.11e-6;%Vm;
Vmeps=Vm;%8.6e-6;%Vm;
Vmeta=Vm;%10.6e-6;%Vm;
Vmsn=Vm;
%% Thermodynamic data from Shim KOLLA HUR DET BLIR MED DU2009
if T<1357.77
    G0alphacu=-7770.458+130.485235*T-24.112392*T*log(T)-0.00265684*T^2+1.29223E-7*T^3+52478*T^-1;
else
    error('Temperature not implemented, G0alphaCu')
end
if T<505.078
    G0bctsn=-5855.135+65.443315*T-15.961*T*log(T)-18.8702e-3*T^2+3.121167E-6*T^3-61960*T^-1; %GHSERSN
elseif T<800
    %     error('Temperature not implemented, G0bctSn')
    G0bctsn=2524.724+4.005269*T-8.2590486*T*log(T)-16.814429e-3*T^2+2.623131e-6*T^3-1081244*T^-1-1.2307e25*T^-9;
end
if T<505.078
    %     G0alphasn=-345.135+56.983315*T-15.961*T*log(T)-0.0188702*T^2+3.121167e-6*T^3-61960*T^-1;%-1705.135+60.243315*T-15.961*T*log(T)-18.8702E-3*T^2+3.121167E-6*T^3-61960*T^-1;%-345.135+56.983315*T-15.961*T*log(T)-0.0188702*T^2+3.121167e-6*T^3-61960*T^-1;
    G0alphasn=-1705.135+60.243315*T-15.961*T*log(T)-18.8702E-3*T^2+3.121167E-6*T^3-61960*T^-1;
elseif T<800
    G0alphasn=6674.724-1.194731*T-8.2590486*T*log(T)-16.814429e-3*T^2+2.623131e-6*T^3-1081244*T^-1-1.2307e25*T^-9;
    %     error('Temperature not implemented, G0alphaSn')
end
% Gculiq=5194.277+120.973331*T-24.112392*T*log(T)-2.65684e-3*T^2+0.129223e-6*T^3+52478*T^-1-5.849e-21*T^7;
% Gsnliq=1247.957+51.355548*T-15.961*T*log(T)-18.8702e-3*T^2+3.121167e-6*T^3-61960*T^-1+1.47031e-18*T^7;
% L0liq=-9002.8-5.8381*T;
% L1liq=-20100.4+3.6366*T;
% L2liq=-10528.4;
% 
% L0alpha=-10672-1.4837*T;%-10309.845+1.157*T;%-10672-1.4837*T;%-10309.845+1.158*T;
% L1alpha=-15331.3+6.9539*T;%-16190.032+6.495*T;%-15331.3+6.9539*T;%-16190.032+6.495*T;

% From Du2009
L0liq=-8266.6-6.9973*T;
L1liq=-21662+8.4655*T;
L2liq=-26957.2+12.8887*T;
L3liq=-13222.7+10.0420*T;
G0alphasn=-345.135+56.983315*T-15.961*T*log(T)-0.0188702*T^2+3.121167e-6*T^3-61960*T^-1; %fcc
L0alpha=-11215.8+1.0978*T;
L1alpha=-12884.1+5.4140*T;
Gbctcu=G0alphacu+5888;
%% Gibbs energies, from Shim
syms c1 c2 cs

%From Du2009
Gscu=(1-cs)*G0alphacu+cs*G0alphasn+R*T*((1-cs)*log(1-cs)+cs*log(cs))+cs*(1-cs)*(L0alpha+L1alpha*(1-2*cs));
mucus=-G0alphacu+G0alphasn+R*T*(log(cs)-log(1-cs))+L0alpha*(1-2*cs)+L1alpha*(6*cs^2-6*cs+1);
d2cu=R*T/(cs-cs^2)-2*L0alpha+L1alpha*(12*cs-6);
Gssn=(1-cs)*Gbctcu+cs*G0bctsn+R*T*((1-cs)*log(1-cs)+cs*log(cs));
musns=-Gbctcu+G0bctsn+R*T*(log(cs)-log(1-cs));
d2sn=R*T/(cs-cs^2);

% %Cu6Sn5 with added c-dependence
etat=2e5;%2e5;
etatc=0.435;
Gseta=etat*(cs-etatc)^2-6869.5-0.1589*T+0.455*G0bctsn+0.545*G0alphacu;
muetas=2*etat*(cs-etatc);
d2eta=2*etat;
% %Cu3Sn with added c-dependence
epst=2e5;
epstc=0.249;
Gseps=epst*(cs-epstc)^2-8194.2-0.2043*T+0.25*G0bctsn+0.75*G0alphacu;
mueps=2*epst*(cs-epstc);
d2eps=2*epst;

%% find equilibrium concentrations
% Cu - eta
init_guess=[0,0.3;0.3,0.5];
e1=subs(Gscu,cs,c1)-subs(mucus,cs,c1)*c1==subs(Gseta,cs,c2)-subs(muetas,cs,c2)*c2;
e2=subs(mucus,cs,c1)==subs(muetas,cs,c2);

[ccu_eta,ceta_cu]=vpasolve([e1,e2],[c1,c2],init_guess);

t0=subs(Gscu,cs,ccu_eta)-subs(mucus,cs,ccu_eta)*ccu_eta+subs(mucus,cs,ccu_eta)*cs;

% Cu - eps
init_guess=[0,0.5;0.,0.4];
e1=subs(Gscu,cs,c1)-subs(mucus,cs,c1)*c1==subs(Gseps,cs,c2)-subs(mueps,cs,c2)*c2;
e2=subs(mucus,cs,c1)==subs(mueps,cs,c2);

[chat_cu,ceps_cu]=vpasolve([e1,e2],[c1,c2],init_guess);

t1=subs(Gscu,cs,chat_cu)-subs(mucus,cs,chat_cu)*chat_cu+subs(mucus,cs,chat_cu)*cs;

% eps-eta
init_guess=[0.2,0.4;0.3,0.5];
e1=subs(Gseps,cs,c1)-subs(mueps,cs,c1)*c1==subs(Gseta,cs,c2)-subs(muetas,cs,c2)*c2;
e2=subs(mueps,cs,c1)==subs(muetas,cs,c2);
[ceps_eta,ceta_eps]=vpasolve([e1,e2],[c1,c2],init_guess);

chat_eps = 0.5*(ceps_cu+ceps_eta);

% cepseta=0.21;
% cetaeps=0.45;
t2=subs(Gseps,cs,ceps_eta)-subs(mueps,cs,ceps_eta)*ceps_eta+subs(mueps,cs,ceps_eta)*cs;

%eta-Sn
init_guess=[0.,1;0.,1];
e1=subs(Gssn,cs,c1)-subs(musns,cs,c1)*c1==subs(Gseta,cs,c2)-subs(muetas,cs,c2)*c2;
e2=subs(musns,cs,c1)==subs(muetas,cs,c2);
[chat_sn,ceta_sn]=vpasolve([e1,e2],[c1,c2],init_guess);
t3=subs(Gssn,cs,chat_sn)-subs(musns,cs,chat_sn)*chat_sn+subs(musns,cs,chat_sn)*cs;

chat_eta=0.5*(ceta_eps+ceta_sn);

%Sn-Cu
init_guess=[0,1;0,1];
e1=subs(Gssn,cs,c1)-subs(musns,cs,c1)*c1==subs(Gscu,cs,c2)-subs(mucus,cs,c2)*c2;
e2=subs(musns,cs,c1)==subs(mucus,cs,c2);
[csn_cu,ccu_sn]=vpasolve([e1,e2],[c1,c2],init_guess);
t4=subs(Gssn,cs,csn_cu)-subs(musns,cs,csn_cu)*csn_cu+subs(musns,cs,csn_cu)*cs;

%Cu3Sn - Sn
init_guess=[0,1;0,1];
e1=subs(Gssn,cs,c1)-subs(musns,cs,c1)*c1==subs(Gseps,cs,c2)-subs(mueps,cs,c2)*c2;
e2=subs(musns,cs,c1)==subs(mueps,cs,c2);
[csn_eps,ceps_sn]=vpasolve([e1,e2],[c1,c2],init_guess);
t5=subs(Gssn,cs,csn_eps)-subs(musns,cs,csn_eps)*csn_eps+subs(musns,cs,csn_eps)*cs;
%% Approximation
%Cu
Agcu=subs(d2cu,cs,chat_cu);
Bgcu=subs(mucus,cs,chat_cu);
Cgcu=subs(Gscu,cs,chat_cu);
Afcu=1/Vmcu*Agcu;
Bfcu=1/Vmcu*Bgcu;
Cfcu=1/Vmcu*Cgcu;
%Sn
Agsn=subs(d2sn,cs,chat_sn);
Bgsn=subs(musns,cs,chat_sn);
Cgsn=subs(Gssn,cs,chat_sn);
Afsn=1/Vmsn*Agsn;
Bfsn=1/Vmsn*Bgsn;
Cfsn=1/Vmsn*Cgsn;
%Cu6Sn5
Ageta=subs(d2eta,cs,chat_eta);
Bgeta=subs(muetas,cs,chat_eta);
Cgeta=subs(Gseta,cs,chat_eta);
Afeta=1/Vmeta*Ageta;
Bfeta=1/Vmeta*Bgeta;
Cfeta=1/Vmeta*Cgeta;

%Cu3Sn
Ageps=subs(d2eps,cs,chat_eta);
Bgeps=subs(mueps,cs,chat_eta);
Cgeps=subs(Gseps,cs,chat_eta);
Afeps=1/Vmeps*Ageps;
Bfeps=1/Vmeps*Bgeps;
Cfeps=1/Vmeps*Cgeps;


% Gibbs energy
gcu=0.5*Agcu*(cs-chat_cu)^2+Bgcu*(cs-chat_cu)+Cgcu;
mucu=Agcu*(cs-chat_cu)+Bgcu;
geps=0.5*Ageps*(cs-chat_eps)^2+Bgeps*(cs-chat_eps)+Cgeps;
mueps=Ageps*(cs-chat_eps)+Bgeps;
geta=0.5*Ageta*(cs-chat_eta)^2+Bgeta*(cs-chat_eta)+Cgeta;
mueta=Ageta*(cs-chat_eta)+Bgeta;
gsn=0.5*Agsn*(cs-chat_sn)^2+Bgsn*(cs-chat_sn)+Cgsn;
musn=Agsn*(cs-chat_sn)+Bgsn;
% Free energy density
fcu=0.5*Afcu*(cs-chat_cu)^2+Bfcu*(cs-chat_cu)+Cfcu;
feps=0.5*Afeps*(cs-chat_eps)^2+Bfeps*(cs-chat_eps)+Cfeps;
feta=0.5*Afeta*(cs-chat_eta)^2+Bfeta*(cs-chat_eta)+Cfeta;
fsn=0.5*Afsn*(cs-chat_sn)^2+Bfsn*(cs-chat_sn)+Cfsn;
mucu = mucu/Vmcu;
mueps = mueps/Vmeps;
mueta = mueta/Vmeta;
musn = musn/Vmsn;

%% plot
%Gibbs
figure()
hold on
p1=fplot(gcu,[0,1]);
% set(p1,'Color','blue')
p2=fplot(geta,[0,1]);
% set(p2,'Color','red');
p3=fplot(gsn,[0,1]);
% set(p3,'Color','green')
p7=fplot(geps,[0,1]);
% set(p2,'Color','magenta');
p4=fplot(Gscu,[0,1]);
% set(p4,'Color','black')
p5=fplot(Gssn,[0,1]);
p6=fplot(Gseta,[0,1]);
% set(p4,'Color','black')
% legend('Cu','eps','eta','Sn','Gcu','Gsn')
% p6=ezplot(t1,[0,1]);
% p7=ezplot(t3,[0,1]);
% p8=ezplot(t4,[0,1]);
% set(p6,'Color','black')
% set(p7,'Color','black')
% set(p8,'Color','black')
title(' ')
ylabel('Gibbs energy [J/mol]')
xlabel('Mole fraction of Sn')
axis([0,1,-4e4,0]);
legend('Cu-fit','Cu6Sn5-fit','Sn-fit','Cu','Sn','Location','sw')
hold off

%Free energy density
figure()
hold on
p1=fplot(fcu,[0,1]);
% set(p1,'Color','blue')
p7=ezplot(feps,[0,1]);
p2=fplot(feta,[0,1]);
% set(p2,'Color','red');
p3=fplot(fsn,[0,1]);
% set(p3,'Color','green')
% set(p2,'Color','magenta');
p4=fplot(1/Vmcu*Gscu,[0,1]);
% set(p4,'Color','black')
p5=fplot(1/Vmsn*Gssn,[0,1]);
% set(p4,'Color','black')
legend('Cu','eps','eta','Sn','Gcu','Gsn')
title('Free energy density')
axis([0,1,-2e9,0]);
hold off


figure()
hold on
p1=fplot(mucu,[0,1]);
% set(p1,'Color','blue')
p4=fplot(mueps,[0,1]);
p2=fplot(mueta,[0,1]);
% set(p2,'Color','red');
p3=fplot(musn,[0,1]);
legend('cu','cu3sn','cu6sn5','sn')
hold off



disp('Cu')
Afcu=eval(Afcu)
Bfcu=eval(Bfcu)
Cfcu=eval(Cfcu)
disp('--------')
disp('Cu3Sn')
Afeps=eval(Afeps)
Bfeps=eval(Bfeps)
Cfeps=eval(Cfeps)
disp('--------')
disp('Cu6Sn5')
Afeta=eval(Afeta)
Bfeta=eval(Bfeta)
Cfeta=eval(Cfeta)
disp('---------')
disp('Sn')
Afsn=eval(Afsn)
Bfsn=eval(Bfsn)
Cfsn=eval(Cfsn)
disp('-----------')
disp('Interface concentrations')
chatcu=chat_cu
chateps = chat_eps
chateta=chat_eta
chatsn=chat_sn

ccu_eps=eval(chat_cu)
ceps_cu = eval(ceps_cu)
ceta_eps=eval(ceta_eps)
ceps_eta=eval(ceps_eta)
csn_eta=eval(chat_sn)
ceta_sn=eval(ceta_sn)

ccu_eta=eval(ccu_eta)
ceta_cu=eval(ceta_cu)
ccu_sn=eval(ccu_sn)
csn_cu=eval(csn_cu)
csn_eps=eval(csn_eps)
ceps_sn=eval(ceps_sn)

