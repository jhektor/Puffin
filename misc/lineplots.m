% fname = 'moelans2011fig2_Lmean_Mvar_sharp.csv'
% fname = 'moelans2011fig3_Lvar_Mvar_sharp_barrier3.csv'

% fname = 'line220C-lennartfig32-34.csv'
% fname = 'line220C-lennartfig32-34-Aeps2e6Aeta2e6-D12.csv'
fname = '3d-1x600x1-25nm-bothIMC-0.02-50h-bdf2.csv'
% fname = '3d-1x180x1-25nm-bothIMC-0.02-50h-v2.csv'
% fname = 'line-50h-bothnuc-D16-fine-lunarc2.csv'

% fname = 'RT_slow-50nm-wi0-350nm-IC1-Limcimccorr-Lbulk0.csv'



dy = 25; %widht of mesh [nm]
dz = 25;

data = csvread(fname,1,0); %read csv file skipping the first row

time = data(:,1)/3600;
cu = data(:,2)/dy/dz/1000; %thickness in um
eps = data(:,3)/dy/dz/1000; 
eta = data(:,4)/dy/dz/1000;
sn = data(:,5)/dy/dz/1000;

stepsize = data(:,6);
energy = data(:,8);

figure()
plot(time,cu,time,eta,time,sn);
legend('Cu','Cu6Sn5','Sn');
ylabel('Thickness [um]')
xlabel('Time [h]')

figure()
plot(time,eps./(eps+eta),time,eta./(eps+eta))
legend('Cu3Sn','Cu6Sn5')


figure()
xdata = [4, 8, 16, 48];
totdata = [1.7, 6.4, 9.8, 13.1];
etadata = [1.7, 6.4, 7.2, 9.4];
epsdata = [0, 0, 2.6, 3.7];
plot(time,(eps+eta),'r',xdata,totdata,'ro',time,eps,'b',xdata,epsdata,'bo',time,eta,'g',xdata,etadata,'go')

legend('sim total','exp total', 'sim Cu3Sn', 'exp Cu3Sn', 'sim Cu6Sn5', 'exp Cu6Sn5','location','best')
ylabel('IMC thickness [um]')
xlabel('Time [h]')



figure()
plot(time,energy)
xlabel('Time [s]')
ylabel('Total energy [eV]')
% 
figure()
plot(time*3600,stepsize)
xlabel('Time [h]')
ylabel('Size of timestep [s]')











% imc = (data(:,3)+data(:,4))/dy; %half thickness of imc layer [micrometer]
% % imc(1)=1.5;
% cu = data(:,2)/dy;%/1000;
% sn = data(:,5)/dy;%/1000;
% figure()
% plot(sqrt(time),cu,sqrt(time),imc,sqrt(time),sn)
% legend('cu','imc','sn')
% figure()
% plot((time),imc)
% energy = data(:,8); %total energy of the system [eV]
% % energy(1) = data(2,8); 
% step_size = data(:,6);
% 
% % fname2 = 'moelans2011fig3_Limc_sn_Msn_sharp.csv'
% % data2 = csvread(fname2,1,0); %read csv file skipping the first row
% % time2 = data2(:,1);
% % imc2 = 0.5*data2(:,3)/dy/1000; %half thickness of imc layer [micrometer]
% % imc2(1)=1.5;
% % fname3 = 'moelans2011fig3_Lcu_imc_Mimc_sharp.csv'
% % data3 = csvread(fname3,1,0); %read csv file skipping the first row
% % time3 = data3(:,1);
% % imc3 = 0.5*data3(:,3)/dy/1000; %half thickness of imc layer [micrometer]
% % imc3(1)=1.5;
% 
%  k=0.5*0.0045;
%  y = k*sqrt(time)+0;
% % 
% % %linear fit to my simulation
% fit_type = fittype({'x'}); %fits data to y = kx
% fit_object = fit(sqrt(time),imc/1000,fit_type);
% figure()
% plot(fit_object,'-k')
% hold on
%  plot(sqrt(time),imc/1000,'r-',sqrt(time),y,'--k');
% % xlabel('Time [s^{1/2}]')
% % ylabel('Half thickness [um]')
% % legend('varying M L_{imc-sn}','M_{sn} L_{imc-sn}','M_{imc} L_{cu-imc}','moelans','location','southeast');
% % axis([0,4000,1.5,10])
% 
% figure()
% plot(time,energy)
% xlabel('Time [s]')
% ylabel('Total energy [eV]')
% 
% figure()
% plot(time,step_size)
% xlabel('Time [s]')
% ylabel('Size of timestep [s]')
% 






