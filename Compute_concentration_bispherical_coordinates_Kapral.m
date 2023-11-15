%compute the concentration field from [Reigh and Kapral, Soft Matter, 2015, 11,
%3149], and Michelin, Lauga : Autophoretic locomotion from geometric asymmetry
%% initialize
clear all;

writedir='Z:\PROJETS\2022_tracers\calculations\Full calculations of hyd fields around dimers\vectorplot_freespace\Mob';% where to write the data
writedata='yes';% save the data in folder writedir?
computestreamF='yes';% compute the stream function or just the concentration field ?

SS=[];
LL=[];
   
AR=[0.54,0.84,1.21,1.55,1.85,2.18,2.50,1]; % aspect ratio between the passive and catalytic spheres.


for AA=1:length(AR)

    % INITIALIZE THE PARAMETERS
    % -------------------------
    lf=10;% minimum 10
    L=[0:lf]';
    blim=15; % limit distance to compute the velocity field and concentration.
    step=0.1;% distance increment

    %0<theta<pi
    %mu=cos(theta) between -1 and 1
    %0<phi<2pi
    r1=1;%radius of sphere 1 (catalytic)
    r2=r1*AR(AA);%radius of sphere 2 (non-catalytic)
    Mobility1=1;%diffusiophoretic mobility of the catalytic part
    Mobility2=1;%diffusiophoretic mobility of the non-catalytic part.


    d=r1+r2+r1/10; %distance between the two spheres
    a=(1/(2*d))*sqrt(((d^2-r1^2-r2^2)^2-4*r1^2*r2^2));
    tau1=log(a/r1+sqrt(1+(a/r1)^2));
    tau2=-log(a/r2+sqrt(1+(a/r2)^2));

    A=-1; % Activity. Negative for repulsive interaction, positive for attractive interaction.


    % COMPUTE THE COEFFICIENTS
    % ------------------------

    M=Mij(tau1,lf,a,A);
    N=Nij(tau1,lf,a,A);
    O=Oij(tau2,lf);
    P=Pij(tau2,lf);
    E=Ei(tau1,a,A,lf);

    % get the Al and Bl coefficients
    Bl=(N-M*(O\P))\E;
    Al=-(O\P)*Bl;


    % find the velocity of the sphere
    L1=[1:lf-1]';
    Speed=compute_speed(L1,tau1,tau2,a,Mobility1,Mobility2,Al,Bl); % compute the phoretic velocity of the dimer
    disp(Speed);

    figure(1)
    clf(1)
    plot(LL,SS,'o');

    % compute the coefficients a,b,c,d
    al=[];bl=[];cl=[];dl=[];

    for l=1:lf-1
            [all,bll,cll,dll]=compute_coeff(l,tau1,tau2,a,Mobility1,Mobility2,Al,Bl,Speed);
            al(l)=all; bl(l)=bll; cl(l)=cll; dl(l)=dll;
    end
    disp([al',bl',cl',dl']);

    %% compute for each position

    Z=[-blim-step/2:step:blim+step/2];
    X=[0-step/2:step:blim+step/2];
    z1=a*coth(tau1);
    z2=a*coth(tau2);
    r1check=a/sinh(tau1);
    r2check=a/sinh(tau2);

    c=zeros(length(Z),length(X));

    TAU=zeros(length(Z),length(X));
    MU=zeros(length(Z),length(X));
    PSY=zeros(length(Z),length(X));% The stream Function

    for tz=1:length(Z)   
        z=Z(tz);
        y=0;
        disp([tz length(Z)]);

        for tx=1:length(X)
            x=X(tx);    
            rho=sqrt(x^2+y^2);


            if((x^2+y^2+(z-z1)^2)>=r1^2 && (x^2+y^2+(z-z2)^2)>=r2^2)
                R=sqrt(z^2+rho^2);
                Q=sqrt((R^2+a^2)^2-(2*a*z)^2);
                tau=asinh(2*a*z/Q);
                sigma=acos((R^2-a^2)/Q);

                if(isequal(x,0))
                    sigma=0;
                end

                mu=cos(sigma);
                MU(tz,tx)=mu;
                TAU(tz,tx)=tau;
                S=(Al.*exp((L+0.5)*(tau))+Bl.*exp(-(L+0.5)*(tau))).*legendreP(L,mu);
                Ssummed=sum(S);

                c(tz,tx)=-(sqrt(cosh(tau)-mu))*(Ssummed);

                if(isequal(computestreamF,'yes'))
                WVsummed=0;

                for l=1:lf-1         %checked
                    Vl=legendrePw(l-1,mu)-legendrePw(l+1,mu);
                    Wl=al(l)*cosh((l-0.5)*tau)+bl(l)*sinh((l-0.5)*tau)+cl(l)*cosh((l+1.5)*tau)+dl(l)*sinh((l+1.5)*tau);
                    WVsummed=WVsummed+Vl*Wl;         
                end

                PSY(tz,tx)=((cosh(tau)-mu)^(-1.5))*WVsummed; % The stream Function
                end

            else
                c(tz,tx)=NaN;
                if(isequal(computestreamF,'yes'))
                    PSY(tz,tx)=NaN;
                end
            end 
        end
    end

    %% ----- plot the results

    if(length(X)==1)
        figure(4)
        clf(4)
        hold on
        plot(abs(Z)',abs(c),'o')
        xlabel( 'z position')
        ylabel('c')
        set(gca(),'XScale','log','YScale','log')
        hold off
        efmenu
    end

    [XX,ZZ]=meshgrid(X,Z);

    figure(5)
    clf(5)
    set(gcf, 'renderer', 'zbuffer')
    hold on
    h =surf(XX,ZZ,real(c));
    colormap(parula)
    shading interp
    xlabel('Z pos (\mu{}m)');
    ylabel('X pos (\mu{}m)');
    title(['concentration Field'])
    cb = colorbar;
    h.FaceColor='flat';
    ylabel(cb, 'concentration (a.u)')
    pbaspect([1 2 1])
    hold off

    figure(6)
    clf(6)
    hold on
    xabs=[0 blim];
    yabs=[-blim blim];
    imagesc(xabs,yabs,real(c));
    pbaspect([1 2 1])
    hold off

    if(isequal(computestreamF,'yes'))
        figure(7)
        clf(7)
        hold on
        h =surf(XX,ZZ,PSY);
        contour(XX,ZZ,PSY);
        colormap(parula)
        shading interp
        xlabel('Z pos (\mu{}m)');
        ylabel('X pos (\mu{}m)');
        title(['Stream function'])
        cb = colorbar;
        h.FaceColor='flat';
        ylabel(cb, 'concentration (a.u)')
        pbaspect([1 2 1])
        hold off

        dPSYx=(PSY(1:end-1,2:end)-PSY(1:end-1,1:end-1))/step;
        dPSYz=(PSY(2:end,1:end-1)-PSY(1:end-1,1:end-1))/step;
        XXv=XX(1:end-1,1:end-1);
        ZZv=ZZ(1:end-1,1:end-1);
        RHO=sqrt(XXv.^2);

        ux=dPSYz./(RHO);% obtain the velocity vector components
        uz=-dPSYx./(RHO);

        startx=[X(1):1:X(end)]';
        startz=[Z(1):1:Z(end)]';
        [SX,SZ]=meshgrid(startx,startz);

        figure(8)
        clf(8)
        hold on
        colormap parula
        q = quiver(XXv,ZZv,ux,uz,1);
        h=streamslice(XXv,ZZv,ux,uz);
        set(h,'LineWidth',2)
        set(h,'Color','r');
        axis([min(X) max(X) min(Z) max(Z)])
        xlabel('Xpos')
        ylabel('Zpos')
        pbaspect([1 2 1])
        hold off
    end


    if(isequal(writedata,'yes'))
        cd(writedir)
        savefig(5,['Fig_concentration_field_' num2str(r2*100) '_Mob' num2str(Mobility1*10)]);   
        savefig(1,['Fig_concentration_plotprofiles_X_Z_' num2str(r2*100) '_Mob' num2str(Mobility1*10)]);   
        saveas(5,['Fig_concentration_field_' num2str(r2*100) '_Mob' num2str(Mobility1*10) '.png']);   
        saveas(1,['Fig_concentration_plotprofiles_X_Z_' num2str(r2*100) '_Mob' num2str(Mobility1*10) '.png']);  

    if(isequal(computestreamF,'yes'))
        savefig(7,['Fig_streamfunction_' num2str(r2*100) '_Mob' num2str(Mobility1*10)]);   
        savefig(8,['Fig_quiverplot_' num2str(r2*100) '_Mob' num2str(Mobility1*10)]);
        saveas(7,['Fig_streamfunction_' num2str(r2*100) '_Mob' num2str(Mobility1*10) '.png']);   
        saveas(8,['Fig_quiverplot_' num2str(r2*100) '_Mob' num2str(Mobility1*10) '.png']);  
        save(['DATA_' num2str(r2*100) '_Mob' num2str(Mobility1*10)],'XX','ZZ','c','XXv','ZZv','ux','uz','PSY','Speed','Mobility1','Mobility2','r1','r2','A','blim','step','lf','d','z1','z2','OMEGA');
    end
    end

end


%we need to compute the matrices M, N, O, P.
%% define functions
function M=Mij(tau11,lf,a,A)
    M=zeros(lf+1,lf+1);
    for n=1:lf+1
        M(n,n)=(sinh(tau11)+(2*n-1)*cosh(tau11))*exp((n-0.5)*tau11);
    end
    
    for n=1:lf
        M(n+1,n)=-n*exp((n-0.5)*tau11);
        M(n,n+1)=-n*exp((n+0.5)*tau11);
    end
end

function N=Nij(tau11,lf,a,A)
    N=zeros(lf+1,lf+1);

    for n=1:lf+1
        N(n,n)=(sinh(tau11)-(2*n-1)*cosh(tau11))*exp(-(n-0.5)*tau11);
    end
    
    for n=1:lf
        N(n+1,n)=n*exp(-(n-0.5)*tau11);
        N(n,n+1)=n*exp(-(n+0.5)*tau11);
    end
end

function O=Oij(tau22,lf)
    O=zeros(lf+1,lf+1);
    for n=1:lf+1
        O(n,n)=(sinh(tau22)+(2*n-1)*cosh(tau22))*exp((n-0.5)*tau22);
    end
    
    for n=1:lf
        O(n+1,n)=-n*exp((n-0.5)*tau22);
        O(n,n+1)=-n*exp((n+0.5)*tau22);
    end
end

function P=Pij(tau22,lf)
    P=zeros(lf+1,lf+1);

    for n=1:lf+1
        P(n,n)=(sinh(tau22)-(2*n-1)*cosh(tau22))*exp(-(n-0.5)*tau22);
    end
    
    for n=1:lf
        P(n+1,n)=n*exp(-(n-0.5)*tau22);
        P(n,n+1)=n*exp(-(n+0.5)*tau22);
    end
end

function E=Ei(tau11,a,A,lf)%checked
    E=zeros(lf+1,1);

    for n=1:lf+1
        E(n)=-2*a*A*sqrt(2)*exp(-(n-0.5)*abs(tau11));
    end
end

function [al,bl,cl,dl]=compute_coeff(l,tau1,tau2,a,M1,M2,AL,BL,V)%checked
    fl=(a^2)*l*(l+1)/(sqrt(2)*(2*l-1)*(2*l+1)*(2*l+3));
    Gl=V*fl;

    S1=cosh((l-0.5)*tau1);
    S2=sinh((l-0.5)*tau1);
    S3=cosh((l+1.5)*tau1);
    S4=sinh((l+1.5)*tau1);
    S=-Gl*((2*l+3)*exp(-(l-0.5)*abs(tau1))-(2*l-1)*exp(-(l+1.5)*abs(tau1)));

    R1=cosh((l-0.5)*tau2);
    R2=sinh((l-0.5)*tau2);
    R3=cosh((l+1.5)*tau2);
    R4=sinh((l+1.5)*tau2);
    R=-Gl*((2*l+3)*exp(-(l-0.5)*abs(tau2))-(2*l-1)*exp(-(l+1.5)*abs(tau2)));

    T1=(2*l-1)*sinh((l-0.5)*tau1);
    T2=(2*l-1)*cosh((l-0.5)*tau1);
    T3=(2*l+3)*sinh((l+1.5)*tau1);
    T4=(2*l+3)*cosh((l+1.5)*tau1);
    T=(2*l-1)*(2*l+3)*Gl*(exp(-(l-0.5)*abs(tau1))-exp(-(l+1.5)*abs(tau1)))+2*a*M1*Phil(l,tau1,AL,BL);

    Q1=(2*l-1)*sinh((l-0.5)*tau2);
    Q2=(2*l-1)*cosh((l-0.5)*tau2);
    Q3=(2*l+3)*sinh((l+1.5)*tau2);
    Q4=(2*l+3)*cosh((l+1.5)*tau2);
    Q=-(2*l-1)*(2*l+3)*Gl*(exp(-(l-0.5)*abs(tau2))-exp(-(l+1.5)*abs(tau2)))+2*a*M2*Phil(l,tau2,AL,BL);

    %we have all the coeff for the set of linear equations. Now we want to
    %solve for al, bl, cl, dl.
    M=[S1 S2 S3 S4;R1 R2 R3 R4;T1 T2 T3 T4;Q1 Q2 Q3 Q4];
    B=[S;R;T;Q];
    Ma=M; Ma(:,1)=B;
    Mb=M; Mb(:,2)=B;
    Mc=M; Mc(:,3)=B;
    Md=M; Md(:,4)=B;

    detM=det(M);
    al=det(Ma)/detM;
    bl=det(Mb)/detM;
    cl=det(Mc)/detM;
    dl=det(Md)/detM;
end

function P=Phil(l,tau,AL,BL)%checked
    Alm1=AL(l);
    Alp1=AL(l+2);
    Al=AL(l+1);

    Blm1=BL(l);
    Blp1=BL(l+2);
    Bl=BL(l+1);

    p1=-l*(l+1)/(2*(2*l+1));
    p2=Alm1*exp((l-0.5)*tau)-Al*(2*cosh(tau))*exp((l+0.5)*tau)+Alp1*exp((l+1.5)*tau)+Blm1*exp(-(l-0.5)*tau)-Bl*(2*cosh(tau))*exp(-(l+0.5)*tau)+Blp1*exp(-(l+1.5)*tau);
    P=p1*p2;
end


function Speed=compute_speed(L,tau1,tau2,a,M1,M2,AL,BL)

    SUM1=0;
    SUM2=0;

    for i=1:length(L)
        l=L(i);
        fl=(a^2)*l*(l+1)/(sqrt(2)*(2*l-1)*(2*l+1)*(2*l+3));

        S1=cosh((l-0.5)*tau1);
        S2=sinh((l-0.5)*tau1);
        S3=cosh((l+1.5)*tau1);
        S4=sinh((l+1.5)*tau1);
        S=-((2*l+3)*exp(-(l-0.5)*abs(tau1))-(2*l-1)*exp(-(l+1.5)*abs(tau1)));

        R1=cosh((l-0.5)*tau2);
        R2=sinh((l-0.5)*tau2);
        R3=cosh((l+1.5)*tau2);
        R4=sinh((l+1.5)*tau2);
        R=-((2*l+3)*exp(-(l-0.5)*abs(tau2))-(2*l-1)*exp(-(l+1.5)*abs(tau2)));

        T1=(2*l-1)*sinh((l-0.5)*tau1);
        T2=(2*l-1)*cosh((l-0.5)*tau1);
        T3=(2*l+3)*sinh((l+1.5)*tau1);
        T4=(2*l+3)*cosh((l+1.5)*tau1);
        Tg=(2*l-1)*(2*l+3)*(exp(-(l-0.5)*abs(tau1))-exp(-(l+1.5)*abs(tau1)));
        Tp=2*a*M1*Phil(l,tau1,AL,BL);

        Q1=(2*l-1)*sinh((l-0.5)*tau2);
        Q2=(2*l-1)*cosh((l-0.5)*tau2);
        Q3=(2*l+3)*sinh((l+1.5)*tau2);
        Q4=(2*l+3)*cosh((l+1.5)*tau2);
        Qg=-(2*l-1)*(2*l+3)*(exp(-(l-0.5)*abs(tau2))-exp(-(l+1.5)*abs(tau2)));
        Qp=2*a*M2*Phil(l,tau2,AL,BL);

        %we have all the coeff for the set of linear equations. Now we want to
        %solve for al, bl, cl, dl.
        M=[S1 S2 S3 S4;R1 R2 R3 R4;T1 T2 T3 T4;Q1 Q2 Q3 Q4];
        detM=det(M);

        % part of DetM*a with G as prefactor
        aG=S*(R2*T3*Q4+R3*T4*Q2+R4*T2*Q3-R2*T4*Q3-R4*T3*Q2-R3*T2*Q4)-S2*(R*T3*Q4+R3*T4*Qg+R4*Tg*Q3-R*T4*Q3-R3*Tg*Q4-R4*T3*Qg)+S3*(R*T2*Q4+R2*T4*Qg+R4*Tg*Q2-R4*T2*Qg-T4*Q2*R-Q4*R2*Tg)-S4*(R*T2*Q3+R2*T3*Qg+R3*Tg*Q2-R*T3*Q2-R2*Tg*Q3-R3*T2*Qg);
        % part of detM*a with the velocity slip component
        aP=-S2*(R3*T4*Qp+R4*Tp*Q3-R3*Tp*Q4-R4*T3*Qp)+S3*(R2*T4*Qp+R4*Tp*Q2-R4*T2*Qp-Q4*R2*Tp)-S4*(R2*T3*Qp+R3*Tp*Q2-R2*Tp*Q3-R3*T2*Qp);

        % part of DetM*c with G as prefactor
        cG=S1*(R2*Tg*Q4+R*T4*Q2+R4*T2*Qg-R2*T4*Qg-R4*Tg*Q2-R*T2*Q4)-S2*(R1*Tg*Q4+R*T4*Q1+R4*T1*Qg-R1*T4*Qg-R*T1*Q4-R4*Tg*Q1)+S*(R1*T2*Q4+R2*T4*Q1+R4*T1*Q2-R4*T2*Q1-T4*Q2*R1-Q4*R2*T1)-S4*(R1*T2*Qg+R2*Tg*Q1+R*T1*Q2-R1*Tg*Q2-R2*T1*Qg-R*T2*Q1);
        % part of detM*c with the velocity slip component
        cP=S1*(R2*Tp*Q4+R4*T2*Qp-R2*T4*Qp-R4*Tp*Q2)-S2*(R1*Tp*Q4+R4*T1*Qp-R1*T4*Qp-R4*Tp*Q1)-S4*(R1*T2*Qp+R2*Tp*Q1-R1*Tp*Q2-R2*T1*Qp);

        SUM1=SUM1+(2*l+1)*(aP+cP)/detM;
        SUM2=SUM2+(2*l+1)*fl*(aG+cG)/detM;
    end

    Speed=-SUM1/SUM2;
end

function P=legendrePw(l,x) % manual entry of the legendre polynomials to go faster
    switch l
        case 0
            P=1;
        case 1
            P=x;
        case 2
            P=0.5*(3*x^2-1);
        case 3
            P=0.5*(5*x^3-3*x);
        case 4
            P=(1/8)*(35*x^4-30*x^2+3);
        case 5
            P=(1/8)*(63*x^5-70*x^3+15*x);
        case 6
            P=(1/16)*(231*x^6-315*x^4+105*x^2-5);
        case 7
            P=(1/16)*(429*x^7-693*x^5+315*x^3-35*x);
        case 8
            P=(1/128)*(6435*x^8-12012*x^6+6930*x^4-1260*x^2+35);
        case 9
            P=(1/128)*(12155*x^9-25740*x^7+18018*x^5-4620*x^3+315*x);
        case 10
            P=(1/256)*(46189*x^10-109395*x^8+90090*x^6-30030*x^4+3465*x^2-63);
        case 11
            P=(1/11)*(-10*(1/128)*(12155*x^9-25740*x^7+18018*x^5-4620*x^3+315*x)+21*x*(1/256)*(46189*x^10-109395*x^8+90090*x^6-30030*x^4+3465*x^2-63));
        otherwise
            P=lengendreP(l,x);
    end
        
end

