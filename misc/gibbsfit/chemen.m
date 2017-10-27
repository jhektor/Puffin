close all
clear all

T=25+273.15;
c=(0.0:1e-6:1)'; %Sn
R=8.3144621; 

Vm=1.0%16.29e-6; %m^3/mole
Vmcu=Vm;%7.124e-6;%7.11e-6;%Vm;
Vmeps=Vm;%8.6e-6;%Vm;
Vmeta=Vm;%10.6e-6;%Vm;
Vmsn=Vm;
%% Shim 1996 %Gierlotka 2007 %Shim et al. 1996
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
Gculiq=5194.277+120.973331*T-24.112392*T*log(T)-2.65684e-3*T^2+0.129223e-6*T^3+52478*T^-1-5.849e-21*T^7;
Gsnliq=1247.957+51.355548*T-15.961*T*log(T)-18.8702e-3*T^2+3.121167e-6*T^3-61960*T^-1+1.47031e-18*T^7;
L0liq=-9002.8-5.8381*T;
L1liq=-20100.4+3.6366*T;
L2liq=-10528.4;

L0alpha=-10672-1.4837*T;%-10309.845+1.157*T;%-10672-1.4837*T;%-10309.845+1.158*T;
L1alpha=-15331.3+6.9539*T;%-16190.032+6.495*T;%-15331.3+6.9539*T;%-16190.032+6.495*T;

%% find equilibrium concentrations
syms c1 c2 cs
Gscu=1/Vmcu*((1-cs)*G0alphacu+cs*G0alphasn+R*T*((1-cs)*log(1-cs)+cs*log(cs))+cs*(1-cs)*(L0alpha+L1alpha*(1-2*cs)));
mucus=1/Vmcu*(-G0alphacu+G0alphasn+R*T*(log(cs)-log(1-cs))+L0alpha*(1-2*cs)+L1alpha*(6*cs^2-6*cs+1));
d2cu=1/Vmcu*(R*T/(cs-cs^2)-2*L0alpha+L1alpha*(12*cs-6));
%Sn in Bct
Gssn=1/Vmsn*(cs*G0bctsn+R*T*((1-cs)*log(1-cs)+cs*log(cs)));
musns=1/Vmsn*(G0bctsn+R*T*(log(cs)-log(1-cs)));
d2sn=1/Vmsn*(R*T/(cs-cs^2));
% Sn in liq
% Gssn=1/Vmsn*((1-cs)*Gculiq+cs*Gsnliq+R*T*((1-cs)*log(1-cs)+cs*log(cs))+cs*(1-cs)*(L0liq+L1liq*(1-2*cs)+L2liq*(1-2*cs)^2));
% musns=1/Vmsn*(-Gculiq+Gsnliq+R*T*(log(cs)-log(1-cs))+L0liq*(1-2*cs)+L1liq*(6*cs^2-6*cs+1)+L2liq*(24*cs^2-16*cs^3-10*cs+1));
% d2sn=1/Vmsn*(R*T/(cs-cs^2)-2*L0alpha+L1alpha*(12*cs-6)+L2liq*(48*cs-48*cs^2-10));
%eta
% Gseta=1/Vmeta*(2E5*(cs-0.435)^2-6869.5-0.1589*T+0.455*G0bctsn+0.545*G0alphacu);
Gseta=1/Vmeta*(2E5*(cs-0.455)^2-7129.7+0.4059*T+0.455*G0bctsn+0.545*G0alphacu);
muetas=1/Vmeta*(4E5*(cs-0.455));
d2eta=4E5/Vmeta;
%eps
Gseps=1/Vmeps*(2E6*(cs-0.249).^2+0.75*G0alphacu+0.25*G0bctsn-8194.2-0.2043*T);
muepss=1/Vmeps*(4E6*(cs-0.249));
d2eps=4E6/Vmeps;
% % Cu-eps
% init_guess=[0,0.5;0,0.5];
% e1=subs(Gscu,cs,c1)-subs(mucus,cs,c1)*c1==subs(Gseps,cs,c2)-subs(muepss,cs,c2)*c2;
% e2=subs(mucus,cs,c1)==subs(muepss,cs,c2);
% 
% [chatcu,ccueps]=vpasolve([e1,e2],[c1,c2],init_guess);
% chatcu=0.0143;
% ccueps=0.2;
% t1=subs(Gscu,cs,chatcu)-subs(mucus,cs,chatcu)*chatcu+subs(mucus,cs,chatcu)*cs;

% % Cu-eta
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
% chatsn=0.999;
% cetasn=0.469;
t3=subs(Gssn,cs,chatsn)-subs(musns,cs,chatsn)*chatsn+subs(musns,cs,chatsn)*cs;

figure
hold on
ezplot(mucus,[0,1]);
ezplot(muetas,[0,1]);
ezplot(musns,[0,1]);
hold off
figure
hold on
p1=ezplot(Gscu,[0,1]);
set(p1,'Color','blue')
p7=ezplot(Gseps,[0,1]);
set(p7,'Color','cyan')
p2=ezplot(Gseta,[0,1]);
set(p2,'Color','red');
p3=ezplot(Gssn,[0,1]);
set(p3,'Color','green')
p4=ezplot(t1,[0,1]);
set(p4,'Color','black')
% p5=ezplot(t2,[0,1]);
% set(p5,'Color','black')
p6=ezplot(t3,[0,1]);
set(p6,'Color','black')
title(' ')
% axis([0,1,-15e8,-6e8])
hold off

%% Parabolic approximations
% chateps=0.5*(ccueps+cepseta);
chateta=0.5*(ccueta+cetasn);%0.5*(cetaeps+cetasn);0.5*(cetaeps+cetasn);%0.5*(ccueta+cetasn);%0.5*(cetaeps+cetasn);
Acu=subs(d2cu,cs,chatcu);
% Aeps=subs(d2eps,cs,chateps);
Aeta=subs(d2eta,cs,chateta);
Asn=subs(d2sn,cs,chatsn);
Bcu=subs(mucus,cs,chatcu);
% Beps=subs(muepss,cs,chateps);
Beta=subs(muetas,cs,chateta);
Bsn=subs(musns,cs,chatsn);
Ccu=subs(Gscu,cs,chatcu);
% Ceps=subs(Gseps,cs,chateps);
Ceta=subs(Gseta,cs,chateta);
Csn=subs(Gssn,cs,chatsn);

fcuf=0.5*Acu*(cs-chatcu)^2+Bcu*(cs-chatcu)+Ccu;
% fepsf=0.5*Aeps*(cs-chateps)^2+Beps*(cs-chateps)+Ceps;
fetaf=0.5*Aeta*(cs-chateta)^2+Beta*(cs-chateta)+Ceta;
fsnf=0.5*Asn*(cs-chatsn)^2+Bsn*(cs-chatsn)+Csn;

figure
hold on
p4=ezplot(Gscu,[0,1]);
set(p4,'Color','black')
p5=ezplot(Gssn,[0,1]);
set(p5,'Color','magenta')
p1=ezplot(fcuf,[0,1]);
set(p1,'Color','blue')
p2=ezplot(fetaf,[0,1]);
set(p2,'Color','red');
p3=ezplot(fsnf,[0,1]);
set(p3,'Color','green')

% p7=ezplot(fepsf,[0,1]);
% set(p7,'Color','yellow')
legend('Cu','Sn','Cu - fit','Cu6Sn5','Sn- fit','Location','SouthEast')
title(' ')
xlabel('Sn concentration')
ylabel('Gibbs energy [J/m^3]')
% axis([0,1,-15e8,-6e8])
hold off
% 
% figure
% plot(c,[Gcu,Geta,Gsn])
% legend('cu','imc','sn','cup','etap','snp')
% xlabel('Concentration of Sn')
% ylabel('Gibbs energy')
% title({'Data from Gierlotka', ['T=',num2str(T-273.15),' C']})
% axis([0,1,-16E8,2E8])


Vm=16.29e-6; %m^3/mole
Vmcu=7.11e-6;%7.11e-6;%Vm;
Vmeps=8.6e-6;%Vm;
Vmeta=10.6e-6;%10.6e-6;%Vm;
Vmsn=Vm;
fcu=1/Vmcu*fcuf;
% feps=1/Vmeps*fepsf;
feta=1/Vmeta*fetaf;
fsn=1/Vmsn*fsnf;
figure
hold on
ezplot(fcu,[0,1]);
% ezplot(feps,[0,1]);
ezplot(feta,[0,1]);
ezplot(fsn,[0,1]);
ezplot(1/Vmcu*Gscu,[0,1]);
ezplot(1/Vmsn*Gssn,[0,1]);
hold off
axis([0,1,-3E9,0.5E9])
legend('Cu','eta','sn','Gcu','Gsn')

disp('Cu')
Acu=eval(Acu)
Bcu=eval(Bcu)
Ccu=eval(Ccu)
disp('--------')
disp('Cu6Sn5')
Aeta=eval(Aeta)
Beta=eval(Beta)
Ceta=Ceta
disp('---------')
disp('Sn')
Asn=eval(Asn)
Bsn=eval(Bsn)
Csn=eval(Csn)
disp('-----------')
disp('Interface concentrations')
chatcu=chatcu
chateta=chateta
chatsn=chatsn
ccueta=eval(ccueta)
cetasn=eval(cetasn)