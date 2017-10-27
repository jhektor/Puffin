%% Fit of Gibbs energy, as in Moelans 2015
close all
clear all

T=25+273.15;%453;%25+273.15;
R=8.3144621;

%% Thermodynamic data from Du2009
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
L0liq=-8266.6-6.9973*T;
L1liq=-21662+8.4655*T;
L2liq=-26957.2+12.8887*T;
L3liq=-13222.7+10.0420*T;
G0alphasn=-345.135+56.983315*T-15.961*T*log(T)-0.0188702*T^2+3.121167e-6*T^3-61960*T^-1; %fcc
L0alpha=-11215.8+1.0978*T;
L1alpha=-12884.1+5.4140*T;
Gbctcu=G0alphacu+5888;



