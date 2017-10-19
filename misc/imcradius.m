fname = 'test2d.csv'

data = csvread(fname,1,0); %read csv file skipping the first row

time = data(:,1);

k=0.04;
y = k*sqrt(time)+0;

imc = sqrt(data(:,2)/pi)/1000; %radius of imc sphere [micrometer]
figure()
plot(sqrt(time),imc,sqrt(time),y)
energy = data(:,5); %total energy of the system [eV]
% energy(1) = data(2,7); 
step_size = data(:,4);




figure()
plot(time,energy)
xlabel('Time [s]')
ylabel('Total energy [eV]')

figure()
plot(time,step_size)
xlabel('Time [s]')
ylabel('Size of timestep [s]')







