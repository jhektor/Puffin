x=0:100:100000;
eta_cu=0.5*(1-tanh(1e-1*(x-48500)));
eta_sn=0.5*(tanh(1e-1*(x-51500)))+0.5;
eta_imc=1-eta_cu-eta_sn;
plot(x,eta_cu,'r-*',x,eta_imc,'g-*',x,eta_sn,'b-*')
hold on