%% Fit of Gibbs energy, as in Moelans 2015
close all
clear all

T=25+273.15;%453;%25+273.15;
R=8.3144621;

Vm=16.29e-6; %m^3/mole
Vmcu=7.124e-6;%7.11e-6;%Vm;
Vmeps=8.6e-6;%Vm;
Vmeta=10.6e-6;%Vm;
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
Gculiq=5194.277+120.973331*T-24.112392*T*log(T)-2.65684e-3*T^2+0.129223e-6*T^3+52478*T^-1-5.849e-21*T^7;
Gsnliq=1247.957+51.355548*T-15.961*T*log(T)-18.8702e-3*T^2+3.121167e-6*T^3-61960*T^-1+1.47031e-18*T^7;
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
% Gscu=((1-cs)*G0alphacu+cs*G0alphasn+R*T*((1-cs)*log(1-cs)+cs*log(cs))+cs*(1-cs)*(L0alpha+L1alpha*(1-2*cs)));
% mucus=(-G0alphacu+G0alphasn+R*T*(log(cs)-log(1-cs))+L0alpha*(1-2*cs)+L1alpha*(6*cs^2-6*cs+1));
% d2cu=(R*T/(cs-cs^2)-2*L0alpha+L1alpha*(12*cs-6));
% %Sn in Bct
% Gssn=(cs*G0bctsn+R*T*((1-cs)*log(1-cs)+cs*log(cs)));
% musns=(G0bctsn+R*T*(log(cs)-log(1-cs)));
% d2sn=(R*T/(cs-cs^2));
% % Sn in liq
% % Gssn=((1-cs)*Gculiq+cs*Gsnliq+R*T*((1-cs)*log(1-cs)+cs*log(cs))+cs*(1-cs)*(L0liq+L1liq*(1-2*cs)+L2liq*(1-2*cs)^2));
% % musns=(-Gculiq+Gsnliq+R*T*(log(cs)-log(1-cs))+L0liq*(1-2*cs)+L1liq*(6*cs^2-6*cs+1)+L2liq*(24*cs^2-16*cs^3-10*cs+1));
% % d2sn=(R*T/(cs-cs^2)-2*L0alpha+L1alpha*(12*cs-6)+L2liq*(48*cs-48*cs^2-10));
%
% %Cu6Sn5
Gseta=-6869.5-0.1589*T+0.455*G0bctsn+0.545*G0alphacu;
% %Cu3Sn
Gseps=-8194.2-0.2043*T+0.25*G0bctsn+0.75*G0alphacu;

%From Du2009
Gscu=(1-cs)*G0alphacu+cs*G0alphasn+R*T*((1-cs)*log(1-cs)+cs*log(cs))+cs*(1-cs)*(L0alpha+L1alpha*(1-2*cs));
mucus=-G0alphacu+G0alphasn+R*T*(log(cs)-log(1-cs))+L0alpha*(1-2*cs)+L1alpha*(6*cs^2-6*cs+1);
d2cu=R*T/(cs-cs^2)-2*L0alpha+L1alpha*(12*cs-6);
Gssn=(1-cs)*Gbctcu+cs*G0bctsn+R*T*((1-cs)*log(1-cs)+cs*log(cs));
musns=-Gbctcu+G0bctsn+R*T*(log(cs)-log(1-cs));
d2sn=R*T/(cs-cs^2);

%% Approximation
% chatcu=0.138074330198467340875139323744; %from old fit
% chatsn=0.9999017777224427579496123643636; %from old fit
% chateta=0.44058110268094034242314439224254; %from old fit
% From Moelans 2015 with Cu3Sn
% chatcu=0.025743;
% chateps=0.25;
% chateta=0.455028;
% chatsn=0.999948;

% From old fits without Cu3Sn
chatcu=0.0059;%25743;%0.1359;%0.13807;
chateps=0.25;
chateta=0.4657;%0.41391;%0.44757;
chatsn=0.99948;%0.99622;%0.999;

res=1;
while res>1e-3
    res
    chatcu0=chatcu;
    chateps0=chateps;
    chateta0=chateta;
    chatsn0=chatsn;
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
    
    %from Moelans2015
    % Acu=2.4e10%1/Vmcu*subs(d2cu,cs,chatcu);
    % Bcu=-4.08e9%1/Vmcu*subs(mucus,cs,chatcu);
    % Ccu=-1.69e9%1/Vmcu*subs(Gscu,cs,chatcu);
    % %Sn
    % Asn=7.19e12%1/Vmsn*subs(d2sn,cs,chatsn);
    % Bsn=4.54e8%1/Vmsn*subs(musns,cs,chatsn);
    % Csn=-2.40e9%1/Vmsn*subs(Gssn,cs,chatsn);
    
    
    %Cu6Sn5 Case 0 in Moelans2015
    Ageta=0.05*Agsn;
    Bgeta=Bgsn;
    Cgeta=Gseta;

    %Check how to do with molar volume
    Afeta=1/Vmeta*Ageta;
    Bfeta=1/Vmeta*Bgeta;
    Cfeta=1/Vmeta*Cgeta;
    %Cu3Sn Case 8
    Ageps=0.05*Agsn;
    Bgeps=Bgcu;
    Cgeps=0.9*Gseps;
    Afeps=1/Vmeps*Ageps;
    Bfeps=1/Vmeps*Bgeps;
    Cfeps=1/Vmeps*Cgeps;
    % Gibbs energy
    gcu=0.5*Agcu*(cs-chatcu)^2+Bgcu*(cs-chatcu)+Cgcu;
    mucu=Agcu*(cs-chatcu)+Bgcu;
    geps=0.5*Ageps*(cs-chateps)^2+Bgeps*(cs-chateps)+Cgeps;
    mueps=Ageps*(cs-chateps)+Bgeps;
    geta=0.5*Ageta*(cs-chateta)^2+Bgeta*(cs-chateta)+Cgeta;
    mueta=Ageta*(cs-chateta)+Bgeta;
    gsn=0.5*Agsn*(cs-chatsn)^2+Bgsn*(cs-chatsn)+Cgsn;
    musn=Agsn*(cs-chatsn)+Bgsn;
    % Free energy density
    fcu=0.5*Afcu*(cs-chatcu)^2+Bfcu*(cs-chatcu)+Cfcu;
    feps=0.5*Afeps*(cs-chateps)^2+Bfeps*(cs-chateps)+Cfeps;
    feta=0.5*Afeta*(cs-chateta)^2+Bfeta*(cs-chateta)+Cfeta;
    fsn=0.5*Afsn*(cs-chatsn)^2+Bfsn*(cs-chatsn)+Cfsn;
    
    
    
    % Find equilibrium concentrations
    %Cu-Cu3Sn
    init_guess=[0,0.5;0.1,0.4];
    e1=subs(gcu,cs,c1)-subs(mucu,cs,c1)*c1==subs(geps,cs,c2)-subs(mueps,cs,c2)*c2;
    e2=subs(mucu,cs,c1)==subs(mueps,cs,c2);
     
    [chatcu,ccueps]=vpasolve([e1,e2],[c1,c2],init_guess);
    t1=subs(gcu,cs,chatcu)-subs(mucu,cs,chatcu)*chatcu+subs(mucu,cs,chatcu)*cs;
    %Cu3Sn-Cu6Sn5
    init_guess=[0.1,0.4;0.3,0.5];
    e1=subs(geps,cs,c1)-subs(mueps,cs,c1)*c1==subs(geta,cs,c2)-subs(mueta,cs,c2)*c2;
    e2=subs(mueps,cs,c1)==subs(mueta,cs,c2);
     
    [cepseta,cetaeps]=vpasolve([e1,e2],[c1,c2],init_guess);
    chateps=ccueps;%0.5*(ccueps+cepseta);
    t2=subs(geps,cs,cepseta)-subs(mueps,cs,cepseta)*cepseta+subs(mueps,cs,cepseta)*cs;
    
    %Cu6Sn5-Sn
    init_guess=[0.3,1;0.3,1];
    e1=subs(gsn,cs,c1)-subs(musn,cs,c1)*c1==subs(geta,cs,c2)-subs(mueta,cs,c2)*c2;
    e2=subs(musn,cs,c1)==subs(mueta,cs,c2);
    [chatsn,cetasn]=vpasolve([e1,e2],[c1,c2],init_guess);
    % chatsn=0.999;
    % cetasn=0.469;
    
    t3=subs(gsn,cs,chatsn)-subs(musn,cs,chatsn)*chatsn+subs(musn,cs,chatsn)*cs;
    chateta=cetasn;%0.5*(cetaeps+cetasn);

%     % Cu-Cu6Sn5
%     init_guess=[0,0.5;0,0.5];
%     e1=subs(gcu,cs,c1)-subs(mucu,cs,c1)*c1==subs(geta,cs,c2)-subs(mueta,cs,c2)*c2;
%     e2=subs(mucu,cs,c1)==subs(mueta,cs,c2);
%     
%     [chatcu,ccueta]=vpasolve([e1,e2],[c1,c2],init_guess);
%     
%     t1=subs(gcu,cs,chatcu)-subs(mucu,cs,chatcu)*chatcu+subs(mucu,cs,chatcu)*cs;
%     
%     %Cu6Sn5-Sn
%     init_guess=[0.3,1;0.3,1];
%     e1=subs(gsn,cs,c1)-subs(musn,cs,c1)*c1==subs(geta,cs,c2)-subs(mueta,cs,c2)*c2;
%     e2=subs(musn,cs,c1)==subs(mueta,cs,c2);
%     [chatsn,cetasn]=vpasolve([e1,e2],[c1,c2],init_guess);
%     % chatsn=0.999;
%     % cetasn=0.469;
%     t3=subs(gsn,cs,chatsn)-subs(musn,cs,chatsn)*chatsn+subs(musn,cs,chatsn)*cs;
%     chateta=0.5*(ccueta+cetasn);

    res=norm([chatcu-chatcu0,chateps-chateps0,chateta-chateta0,chatsn-chatsn0]);
    
end
% geta=0.5*0.01*Agsn*(cs-chateta)^2+Bgeta*(cs-chateta)+Cgeta;
%% plot
%Gibbs
figure()
hold on
p1=ezplot(gcu,[0,1]);
% set(p1,'Color','blue')
ezplot(geps,[0,1]);
p2=ezplot(geta,[0,1]);
% set(p2,'Color','red');
p3=ezplot(gsn,[0,1]);
% set(p3,'Color','green')
% p2=ezplot(geps,[0,1]);
% set(p2,'Color','magenta');
p4=ezplot(Gscu,[0,1]);
% set(p4,'Color','black')
p5=ezplot(Gssn,[0,1]);
% set(p4,'Color','black')
% legend('Cu','eps','eta','Sn','Gcu','Gsn')
p6=ezplot(t1,[0,1]);
p7=ezplot(t3,[0,1]);
p8=ezplot(t2,[0,1]);
set(p6,'Color','black')
set(p7,'Color','black')
set(p8,'Color','black')
title('Gibbs energy')
axis([0,1,-20e3,-5e3]);
hold off

%Free energy density
figure()
hold on
p1=ezplot(fcu,[0,1]);
ezplot(feps,[0,1]);
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

disp('Cu')
Agcu=eval(Agcu)
Bgcu=eval(Bgcu)
Cgcu=eval(Cgcu)
disp('--------')
disp('Cu6Sn5')
Ageta=eval(Ageta)
Bgeta=eval(Bgeta)
Cgeta=Cgeta
disp('---------')
disp('Sn')
Agsn=eval(Agsn)
Bgsn=eval(Bgsn)
Cgsn=eval(Cgsn)
disp('-----------')
disp('Interface concentrations')
chatcu=chatcu
chateta=chateta
chatsn=chatsn
% ccueta=eval(ccueta)
ccueps=eval(ccueps)
cetasn=eval(cetasn)


