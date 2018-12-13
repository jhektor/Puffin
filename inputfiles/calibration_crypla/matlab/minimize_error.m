function err = minimize_error(x,plotflag,lb,ub)
err=0;

E=x(1);
v=x(2);
G=E/(2*(1+v));   %18e9%6e9;%20e9;
K=E/(3*(1-2*v));

mp=[x,0];
mp(1)=G;
mp(2)=K;


lcases=[6];%,2,3];%,2,3];%,3];
emax=0.12;
for i=1:numel(lcases)
    lcase=lcases(i);
    if lcase==1
        S=load('../data/darbandi2012.mat','xx','yy');
%         expstrain=S.xx(2:end);
%         expstress=S.yy(2:end);
        expstrain=S.xx(S.xx<=emax);
        expstress=S.yy(S.xx<=emax);
    elseif lcase==2 %philippi 110
        S=load('../data/philippi2016extra110.mat','xx','yy');
        expstrain=S.xx(S.xx<=emax);
        expstress=S.yy(S.xx<=emax);
    elseif lcase==3 %philippi 001
        S=load('../data/philippi2016extra001.mat','xx','yy');
        expstrain=S.xx(S.xx<=emax);
        expstress=S.yy(S.xx<=emax);
    elseif lcase==4 %Ekinci 110 fig1
        S=load('../data/ekinci2003fig1-110.mat','xx','yy');
        expstrain=S.xx(S.xx<=emax);
        expstress=1e6*S.yy(S.xx<=emax);
    elseif lcase==5 %Kiener2006 Cu
        S=load('../data/kiener2006fig4b.mat','xx','yy');
        expstrain=S.xx(S.xx<=emax);
        expstress=S.yy(S.xx<=emax);
    elseif lcase==6 %Champion2003 Cu polycrystal
        S=load('../data/champion2003fig2cutest.mat','xx','yy');
        expstrain=S.xx(S.xx<=emax);
        expstress=S.yy(S.xx<=emax);
    end
    mp(10)=lcase;
    % write mp to text file so that Fortran kan read it
    save('mp_opt.txt','mp','-ascii','-double')
    % run the simulation in Fortran
    [stat,res]=system('./calib.x','-echo');
    
    % read the results from the mat-file written by Fortran
    load('calib.mat')
    
    mstrain=zeros(size(expstrain));
    mstress=zeros(size(expstress));
    %find simulated strains closest to experimental
    for i=1:size(expstrain,1)
        val=expstrain(i);
        tmp = abs(strain-val);
        [~, idx] = min(tmp);
        mstrain(i)=strain(idx);
        mstress(i)=stress(idx);
    end
    
    %penalty to unrealistic values
    %x0=[G,K,G0,gamma0,mint,hab,QQ,BB,gg];
    
    penalty=1;
    if max(stress)<0
        penalty=1e5
    end
    for i=1:9
        if (x(i)<=lb(i) || x(i)>ub(i)) 
            penalty = 1e5;
        else
            if x(i)>0.5*(lb(i)+ub(i))
                xp = x(i)/ub(i);
                yp = (1 + xp/(1+(10000-1)*(1-xp)))^5;
                penalty = yp*penalty;
                if penalty<0
                    error('penalty negative')
                end
            else
                
                xp = lb(i)/x(i);
                yp = (1 + xp/(1+(10000-1)*(1-xp)))^5;
                penalty = yp*penalty;
                if penalty<0
                    error('penalty negative')
                end
            end
        end
    end
    
    npoints=81;
    weight=1/numel(mstress);
%                  penalty=1;
    err=err+weight*penalty*norm((mstress-expstress)/expstress);
    
    if isnan(err)
        err = 100000;
        'error is nan'
    end
    if plotflag
        colors=['r','g','b','k','c','y'];
        figure(2);
        if lcase==lcases(1)
            clf(2)
        end
        plot(mstrain,mstress,strcat(colors(find(lcases==lcase)),'-'))
        hold on
        plot(expstrain,expstress,strcat(colors(find(lcases==lcase)),'*'))
    end
end

%save('error_and_mps.txt','mp','err','-ascii')
%save('error_and_mps.mat','mp','err','-append')
