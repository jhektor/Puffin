%% Fit of Gibbs energy, as in Moelans 2015
close all
clear all

T=150+273.15;%453;%25+273.15;
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
etat=2e5;
etatc=0.435;
Gseta=etat*(cs-etatc)^2-6869.5-0.1589*T+0.455*G0bctsn+0.545*G0alphacu;
muetas=2*etat*(cs-etatc);
d2eta=2*etat;
% %Cu3Sn with added c-dependence
epst=2e5;
epstc=0.25;
Gseps=epst*(cs-epstc)^2-8194.2-0.2043*T+0.25*G0bctsn+0.75*G0alphacu;
mepss=2*epst*(cs-epstc);
d2eps=2*epst;

%% find equilibrium concentrations
init_guess=[0,0.5;0,0.5];
e1=subs(Gscu,cs,c1)-subs(mucus,cs,c1)*c1==subs(Gseta,cs,c2)-subs(muetas,cs,c2)*c2;
e2=subs(mucus,cs,c1)==subs(muetas,cs,c2);

[chatcu,ccueta]=vpasolve([e1,e2],[c1,c2],init_guess);

t1=subs(Gscu,cs,chatcu)-subs(mucus,cs,chatcu)*chatcu+subs(mucus,cs,chatcu)*cs;

% % eps-eta
% init_guess=[0,0.5;0,0.7];
% e1=subs(Gseps,cs,c1)-subs(muepss,cs,c1)*c1==subs(Gseta,cs,c2)-subs(muetas,cs,c2)*c2;
% e2=subs(muepss,cs,c1)==subs(muetas,cs,c2);
% [cepseta,cetaeps]=vpasolve([e1,e2],[c1,c2],init_guess);
% % cepseta=0.21;
% % cetaeps=0.45;
% t2=subs(Gseps,cs,cepseta)-subs(muepss,cs,cepseta)*cepseta+subs(muepss,cs,cepseta)*cs;

%eta-Sn
init_guess=[0.3,1;0.3,1];
e1=subs(Gssn,cs,c1)-subs(musns,cs,c1)*c1==subs(Gseta,cs,c2)-subs(muetas,cs,c2)*c2;
e2=subs(musns,cs,c1)==subs(muetas,cs,c2);
[chatsn,cetasn]=vpasolve([e1,e2],[c1,c2],init_guess);
t3=subs(Gssn,cs,chatsn)-subs(musns,cs,chatsn)*chatsn+subs(musns,cs,chatsn)*cs;

chateta=0.5*(ccueta+cetasn);

%Sn-Cu
init_guess=[0,1;0,1];
e1=subs(Gssn,cs,c1)-subs(musns,cs,c1)*c1==subs(Gscu,cs,c2)-subs(mucus,cs,c2)*c2;
e2=subs(musns,cs,c1)==subs(mucus,cs,c2);
[csncu,ccusn]=vpasolve([e1,e2],[c1,c2],init_guess);
t4=subs(Gssn,cs,csncu)-subs(musns,cs,csncu)*csncu+subs(musns,cs,csncu)*cs;
%% Approximation
%Cu
Agcu=subs(d2cu,cs,chatcu);
Bgcu=subs(mucus,cs,chatcu);
Cgcu=subs(Gscu,cs,chatcu);
Afcu=1/Vmcu*Agcu;
Bfcu=1/Vmcu*Bgcu;
Cfcu=1/Vmcu*Cgcu;
%Sn
Agsn=subs(d2sn,cs,chatsn);
Bgsn=subs(musns,cs,chatsn);
Cgsn=subs(Gssn,cs,chatsn);
Afsn=1/Vmsn*Agsn;
Bfsn=1/Vmsn*Bgsn;
Cfsn=1/Vmsn*Cgsn;
%Cu6Sn5
Ageta=subs(d2eta,cs,chateta);
Bgeta=subs(muetas,cs,chateta);
Cgeta=subs(Gseta,cs,chateta);
Afeta=1/Vmeta*Ageta;
Bfeta=1/Vmeta*Bgeta;
Cfeta=1/Vmeta*Cgeta;

% Gibbs energy
gcu=0.5*Agcu*(cs-chatcu)^2+Bgcu*(cs-chatcu)+Cgcu;
mucu=Agcu*(cs-chatcu)+Bgcu;
% geps=0.5*Ageps*(cs-chateps)^2+Bgeps*(cs-chateps)+Cgeps;
geta=0.5*Ageta*(cs-chateta)^2+Bgeta*(cs-chateta)+Cgeta;
mueta=Ageta*(cs-chateta)+Bgeta;
gsn=0.5*Agsn*(cs-chatsn)^2+Bgsn*(cs-chatsn)+Cgsn;
musn=Agsn*(cs-chatsn)+Bgsn;
% Free energy density
fcu=0.5*Afcu*(cs-chatcu)^2+Bfcu*(cs-chatcu)+Cfcu;
% feps=0.5*Afeps*(cs-chateps)^2+Bfeps*(cs-chateps)+Cfeps;
feta=0.5*Afeta*(cs-chateta)^2+Bfeta*(cs-chateta)+Cfeta;
fsn=0.5*Afsn*(cs-chatsn)^2+Bfsn*(cs-chatsn)+Cfsn;


%% plot
%Gibbs
figure()
hold on
p1=ezplot(gcu,[0,1]);
% set(p1,'Color','blue')
p2=ezplot(geta,[0,1]);
% set(p2,'Color','red');
p3=ezplot(gsn,[0,1]);
% set(p3,'Color','green')
% p2=ezplot(geps,[0,1]);
% set(p2,'Color','magenta');
p4=ezplot(Gscu,[0,1]);
% set(p4,'Color','black')
p5=ezplot(Gssn,[0,1]);
p6=ezplot(Gseta,[0,1]);
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
axis([0,1,-20e3,-5e3]);
legend('Cu-fit','Cu6Sn5-fit','Sn-fit','Cu','Sn','Location','sw')
hold off

%Free energy density
figure()
hold on
p1=ezplot(fcu,[0,1]);
% set(p1,'Color','blue')
p2=ezplot(feta,[0,1]);
% set(p2,'Color','red');
p3=ezplot(fsn,[0,1]);
% set(p3,'Color','green')
% p2=ezplot(feps,[0,1]);
% set(p2,'Color','magenta');
p4=ezplot(1/Vmcu*Gscu,[0,1]);
% set(p4,'Color','black')
p5=ezplot(1/Vmsn*Gssn,[0,1]);
% set(p4,'Color','black')
% legend('Cu','eps','eta','Sn','Gcu','Gsn')
title('Free energy density')
axis([0,1,-22e8,-5e8]);
hold off


figure()
hold on
p1=ezplot(mucu,[0,1]);
% set(p1,'Color','blue')
p2=ezplot(mueta,[0,1]);
% set(p2,'Color','red');
p3=ezplot(musn,[0,1]);
legend('cu','imc','sn')
hold off



disp('Cu')
Agcu=eval(Afcu)
Bgcu=eval(Bfcu)
Cgcu=eval(Cfcu)
disp('--------')
disp('Cu6Sn5')
Ageta=eval(Afeta)
Bgeta=eval(Bfeta)
Cgeta=Cfeta
disp('---------')
disp('Sn')
Agsn=eval(Afsn)
Bgsn=eval(Bfsn)
Cgsn=eval(Cfsn)
disp('-----------')
disp('Interface concentrations')
chatcu=chatcu
chateta=chateta
chatsn=chatsn
ccueta=eval(ccueta)
cetasn=eval(cetasn)


