%Plots a superellipse
x=-0.5:0.001:0.5;

a=0.01;
b=0.01;
n=2 %0<n<1: concave star; n=1: rhombus; 1<n<2: bulged rhombus; n=2: ellipse; n>2: rounded rectangle

y=b*(1-(x./a).^(n)).^(1/n);


y = 1-0.5*(1+tanh(x/0.1));

plot(x,y,'r',x,-y,'r',-x,y,'r',-x,-y,'r')
axis equal


