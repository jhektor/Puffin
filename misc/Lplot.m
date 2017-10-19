eta=0:0.01:1;
[e1,e2]=meshgrid(eta,eta);
e3=0.5;

mask=e1+e2+e3<=1;

e1=mask.*e1;
e2=mask.*e2;

L12=0.35;
L23=0.705;
L13=0.5;
L=(L12*e1.^2.*e2.^2+L23*e2.^2.*e3.^2+L13*e1.^2.*e3.^2)./(e1.^2.*e2.^2+e2.^2.*e3.^2+e1.^2.*e3.^2);

surf(e1,e2,L,'edgecolor','none')