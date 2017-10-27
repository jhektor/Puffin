close all
clear all

T=523;%125+273.15;
c=(0.0:1e-6:1)'; %Sn
R=8.3144621; 

Vm=16.29e-6;%16.29e-6; %m^3/mole
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
% Gssn=1/Vmsn*(cs*G0bctsn+R*T*((1-cs)*log(1-cs)+cs*log(cs)));
% musns=1/Vmsn*(G0bctsn+R*T*(log(cs)-log(1-cs)));
% d2sn=1/Vmsn*(R*T/(cs-cs^2));
% Sn in liq
Gssn=1/Vmsn*((1-cs)*Gculiq+cs*Gsnliq+R*T*((1-cs)*log(1-cs)+cs*log(cs))+cs*(1-cs)*(L0liq+L1liq*(1-2*cs)+L2liq*(1-2*cs)^2));
musns=1/Vmsn*(-Gculiq+Gsnliq+R*T*(log(cs)-log(1-cs))+L0liq*(1-2*cs)+L1liq*(6*cs^2-6*cs+1)+L2liq*(24*cs^2-16*cs^3-10*cs+1));
d2sn=1/Vmsn*(R*T/(cs-cs^2)-2*L0alpha+L1alpha*(12*cs-6)+L2liq*(48*cs-48*cs^2-10));
%eta
% Gseta=1/Vmeta*(2E5*(cs-0.435)^2-6869.5-0.1589*T+0.455*G0bctsn+0.545*G0alphacu);
Gseta=1/Vmeta*(2E5*(cs-0.455)^2-7129.7+0.4059*T+0.455*G0bctsn+0.545*G0alphacu);
muetas=1/Vmeta*(4E5*(cs-0.455));
d2eta=4E5/Vmeta;
%eps
Gseps=1/Vmeps*(2E6*(cs-0.249).^2+0.75*G0alphacu+0.25*G0bctsn-8194.2-0.2043*T);
muepss=1/Vmeps*(4E6*(cs-0.249));
d2eps=4E6/Vmeps;


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

hold off
% axis([0,1,-2e9,-1.5e9])