T=25+273.15;
R=8.3144598;

%eq_composition of Sn
c_cu = 0.10569; %sn in cu
c_sn = 0.99941; %sn in sn
%Sn in Cu
D0_sn_cu = 2.95e-5;
Q_sn_cu = 177e3;
%Cu in Cu
D0_cu_cu = 3.4e-5;
Q_cu_cu = 195.6e3;

%Sn in Sn
D0_sn_sn = 1.2e-9;
Q_sn_sn = 43.89e3;
%Cu in sn
D0_cu_sn = 2.4e-7;
Q_cu_sn = 33.02e3;

%Interdiffusion coefficients
Dsn_cu=D0_sn_cu*exp(-Q_sn_cu/(R*T));
Dcu_cu=D0_cu_cu*exp(-Q_cu_cu/(R*T));
Dcu_sn=D0_cu_sn*exp(-Q_cu_sn/(R*T));
Dsn_sn=D0_sn_sn*exp(-Q_sn_sn/(R*T));

% Dcu=((1-c_sn)*Dcu_cu+c_cu*Dsn_cu)/((1-c_cu)+c_sn) %????????????
% 
%Cu3Sn
D0_eps = 5.48e-9;
Q_eps = 61.9e3;

%Cu6Sn5
D0_eta = 1.84e-9;
Q_eta = 53.9e3;

Dcu = D0_sn_cu*exp(-Q_sn_cu/(R*T))
Deps=D0_eps*exp(-Q_eps/(R*T))
Deta=D0_eta*exp(-Q_eta/(R*T))
Dsn=D0_sn_sn*exp(-Q_sn_sn/(R*T))