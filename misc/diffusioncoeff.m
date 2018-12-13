D0cusn=2.95e-5;        %Sn in Cu
D0cucu=3.4e-5;         %Cu in Cu
D0eps=5.48e-9;         %Cu3Sn
D0eta=1.84e-9;         %Cu6Sn5
D0snsn=1.2e-9;         %Sn in Sn
D0sncu=2.4e-7;         %Cu in Sn
Qcusn=177.0e3;         %Activation energy [J/mol] Sn in Cu
Qcucu=195.6e3;         %Cu in Cu
Qeps=61.86e3;
Qeta=53.92e3;
Qsnsn=43.89e3;         %Sn in Sn
Qsncu=33.02e3;         %Cu in Sn  

R = 8.314;

T=150+293.15;

Dcu=D0cusn*exp(-Qcusn/(R*T))
Deps=D0eps*exp(-Qeps/(R*T))
Deta=D0eta*exp(-Qeta/(R*T))
Dsn=D0snsn*exp(-Qsnsn/(R*T))