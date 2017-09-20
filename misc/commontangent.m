function f=commontangent(x,ABC1,ABC2,xhat1,xhat2)
A1=ABC1(1);
B1=ABC1(2);
C1=ABC1(3);
A2=ABC2(1);
B2=ABC2(2);
C2=ABC2(3);

f1=0.5*A1*(x(1)-xhat1)^2+B1*(x(1)-xhat1)+C1;
f2=0.5*A2*(x(2)-xhat2)^2+B2*(x(2)-xhat2)+C2;

df1=A1*(x(1)-xhat1)+B1;
df2=A2*(x(2)-xhat2)+B2;

%Residuals
f=[df1-df2;
    f1-df1*x(1)-f2+df2*x(2)];
end