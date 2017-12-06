%Plots a superellipse
x=0:1:2000;
a=500;
b=400;
n=5; %0<n<1: concave star; n=1: rhombus; 1<n<2: bulged rhombus; n=2: ellipse; n>2: rounded rectangle

y=b*(1-(x./a).^(n)).^(1/n);
plot(x,y,'r',x,-y,'r',-x,y,'r',-x,-y,'r')
axis equal


