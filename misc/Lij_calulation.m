syms L12 L23 L13 eta

% Determine L 
% Need to constrain the etas so we can have a single variable
eta1=eta;
eta2=0.5*(1-eta);
eta3=0.5*(1-eta);

f = L12*eta1^2*eta2^2+L23*eta2^2*eta3^2+L13*eta1^2*eta3^2;
g = eta1^2*eta2^2+eta2^2*eta3^2+eta1^2*eta3^2;
L=f/g;
%Now we use L'Hospital to determine the limit of L as eta->1
dfdeta = diff(f,eta)
dgdeta = diff(g,eta)

Ls=simplify(dfdeta/dgdeta)
subs(Ls,eta,1) %Limit as eta -> 1

L=subs(L,L12,0.0711) ;
L=subs(L,L13,0.00119) ;
L=subs(L,L23,0.03) ;

figure()
fplot(L,[-0.1,1.3])
hold on

pf=1e5;
ep=0.01;

L2=(L12*(pf*eta1+ep)^2*(pf*eta2+ep)^2+L13*(pf*eta1+ep)^2*(pf*eta3+ep)^2+L23*(pf*eta2+ep)^2*(pf*eta3+ep)^2)/((pf*eta1+ep)^2*(pf*eta2+ep)^2+(pf*eta1+ep)^2*(pf*eta3+ep)^2+(pf*eta2+ep)^2*(pf*eta3+ep)^2)
L2=subs(L2,L12,0.0711) ;
L2=subs(L2,L13,0.00119) ;
L2=subs(L2,L23,0.03) ;
fplot(L2,[-0.1,1.3])

hold off
figure()
fplot(abs(L-L2),[-0.1,1.3])


