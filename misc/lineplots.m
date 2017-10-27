% fname = 'moelans2011fig2_Lmean_Mvar_sharp.csv'
% fname = 'moelans2011fig3_Lvar_Mvar_sharp_barrier3.csv'
fname = 'test_moelans.csv'


dy = 100; %element height in [nm]


data = csvread(fname,1,0); %read csv file skipping the first row

time = data(:,1);

imc = 0.5*(data(:,3)+data(:,4))/dy; %half thickness of imc layer [micrometer]
% imc(1)=1.5;
cu = data(:,2)/dy;%/1000;
sn = data(:,5)/dy;%/1000;
figure()
plot(sqrt(time),cu,sqrt(time),imc,sqrt(time),sn)
legend('cu','imc','sn')
figure()
plot(sqrt(time),imc)
energy = data(:,8); %total energy of the system [eV]
% energy(1) = data(2,8); 
step_size = data(:,6);

% fname2 = 'moelans2011fig3_Limc_sn_Msn_sharp.csv'
% data2 = csvread(fname2,1,0); %read csv file skipping the first row
% time2 = data2(:,1);
% imc2 = 0.5*data2(:,3)/dy/1000; %half thickness of imc layer [micrometer]
% imc2(1)=1.5;
% fname3 = 'moelans2011fig3_Lcu_imc_Mimc_sharp.csv'
% data3 = csvread(fname3,1,0); %read csv file skipping the first row
% time3 = data3(:,1);
% imc3 = 0.5*data3(:,3)/dy/1000; %half thickness of imc layer [micrometer]
% imc3(1)=1.5;

 k=0.5*0.0045;
 y = k*sqrt(time)+0;
% 
% %linear fit to my simulation
fit_type = fittype({'x'}); %fits data to y = kx
fit_object = fit(sqrt(time),imc/1000,fit_type);
figure()
plot(fit_object,'-k')
hold on
 plot(sqrt(time),imc/1000,'r-',sqrt(time),y,'--k');
% xlabel('Time [s^{1/2}]')
% ylabel('Half thickness [um]')
% legend('varying M L_{imc-sn}','M_{sn} L_{imc-sn}','M_{imc} L_{cu-imc}','moelans','location','southeast');
% axis([0,4000,1.5,10])

figure()
plot(time,energy)
xlabel('Time [s]')
ylabel('Total energy [eV]')

figure()
plot(time,step_size)
xlabel('Time [s]')
ylabel('Size of timestep [s]')







