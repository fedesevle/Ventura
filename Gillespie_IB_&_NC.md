%% Simulador de datos estocásticos IB y NC
clear all;close all
load param.mat;
n=100; %nº de sets
V=1000; % Volumen y Nro de partículas, así la concentración es invariante
resL=40;LL=logspace(0,11,resL); % LL/V tiene que ser misma concentracion que en determinista para comparar
rest=1000; %tiene que ser misma que en determinista para comparar
%% IB
R00j=zeros(rest,resL);RL0j=zeros(rest,resL);R0Lj=zeros(rest,resL);RLLj=zeros(rest,resL);tj=logspace(-4,3,rest);
for m=1:n
    k10=param(m,1);l10=param(m,2);k01=param(m,3);l01=param(m,4);
    for j=1:resL
    R00(1)=V;RL0(1)=0;R0L(1)=0;RLL(1)=0;t(1)=0;L=LL(j)/V;i=0;t(1)=0;t(2)=0; %condiciones iniciales
    tita_final=1/2*(L/(l01/k01+L)+L/(l10/k10+L));contador=0; %condiciones finales
        while contador<V && t(i+1)<max(tj)
            i=i+1;
            dR(i)=(k10*(R00(i)+R0L(i))+k01*(R00(i)+RL0(i)))*L/V+((l10+l01)*RLL(i)+l10*RL0(i)+l01*R0L(i))/V;
            dt=exprnd(1/dR(i));
            t(i+1)=t(i)+dt;
            Pk10_00=k10*R00(i)*L/V/dR(i);
            Pk01_00=k01*R00(i)*L/V/dR(i);
            Pk10_0L=k10*R0L(i)*L/V/dR(i);
            Pk01_L0=k01*RL0(i)*L/V/dR(i);
            Pl10_00=l10*RL0(i)/dR(i);
            Pl01_00=l01*R0L(i)/dR(i);
            Pl10_LL=l10*RLL(i)/dR(i);
            Pl01_LL=l01*RLL(i)/dR(i);
            P=rand();
            if P<Pk10_00
                R00(i+1)=R00(i)-1;
                RL0(i+1)=RL0(i)+1;
                R0L(i+1)=R0L(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00 && P<Pk10_00+Pk01_00
                R00(i+1)=R00(i)-1;
                R0L(i+1)=R0L(i)+1;
                RL0(i+1)=RL0(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00+Pk01_00 && P<Pk10_00+Pk01_00+Pk10_0L
                R0L(i+1)=R0L(i)-1;
                RLL(i+1)=RLL(i)+1;
                RL0(i+1)=RL0(i);
                R00(i+1)=R00(i);                
            elseif P>Pk10_00+Pk01_00+Pk10_0L && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0
                RL0(i+1)=RL0(i)-1;
                RLL(i+1)=RLL(i)+1;
                R0L(i+1)=R0L(i);
                R00(i+1)=R00(i);                
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0 && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00
                RL0(i+1)=RL0(i)-1;
                R00(i+1)=R00(i)+1;
                R0L(i+1)=R0L(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00 && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00
                R0L(i+1)=R0L(i)-1;
                R00(i+1)=R00(i)+1;
                RL0(i+1)=RL0(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00 && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00+Pl10_LL
                RLL(i+1)=RLL(i)-1;
                R0L(i+1)=R0L(i)+1;
                RL0(i+1)=RL0(i);
                R00(i+1)=R00(i);
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00+Pl10_LL
                RLL(i+1)=RLL(i)-1;
                RL0(i+1)=RL0(i)+1;
                R0L(i+1)=R0L(i);
                R00(i+1)=R00(i);
            end
            tita=(RL0(i)+R0L(i)+2*RLL(i))/(2*V);
            if tita_final-tita<0
                contador=contador+1;
            end
        end
        R00j(:,j)=interp1(t,R00,tj);R0Lj(:,j)=interp1(t,R0L,tj);RL0j(:,j)=interp1(t,RL0,tj);RLLj(:,j)=interp1(t,RLL,tj);
        clear R00 RL0 R0L RLL t
    end
    m
    j
    IB(m)=struct('RLL',RLLj,'tita',(2*RLLj+RL0j+R0Lj)/(2*V),'t',tj);
end
%% NC
R00j=zeros(rest,resL);RL0j=zeros(rest,resL);R0Lj=zeros(rest,resL);RLLj=zeros(rest,resL);tj=logspace(-4,3,rest);
for m=1:n
    k=param(m,5);l=param(m,6);w=param(m,7);
    for j=1:resL
    R00(1)=V;RL0(1)=0;R0L(1)=0;RLL(1)=0;t(1)=0;L=LL(j)/V;i=0;t(1)=0;t(2)=0; % condiciones iniciales
    tita_final=(L*l/k+w*L^2)/((l/k)^2+2*L*l/k+w*L^2);contador=0; %condiciones finales
        while contador<V && t(i+1)<max(tj)
            i=i+1;
            dR(i)=(2*k*R00(i)+k*w*(R0L(i)+RL0(i)))*L/V+(2*l*RLL(i)+l*(RL0(i)+R0L(i)))/V;
            dt=exprnd(1/dR(i));
            t(i+1)=t(i)+dt;
            Pk10_00=k*R00(i)*L/V/dR(i);
            Pk01_00=k*R00(i)*L/V/dR(i);
            Pk10_0L=k*w*R0L(i)*L/V/dR(i);
            Pk01_L0=k*w*RL0(i)*L/V/dR(i);
            Pl10_00=l*RL0(i)/dR(i);
            Pl01_00=l*R0L(i)/dR(i);
            Pl10_LL=l*RLL(i)/dR(i);
            Pl01_LL=l*RLL(i)/dR(i);
            P=rand();
            if P<Pk10_00
                R00(i+1)=R00(i)-1;
                RL0(i+1)=RL0(i)+1;
                R0L(i+1)=R0L(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00 && P<Pk10_00+Pk01_00
                R00(i+1)=R00(i)-1;
                R0L(i+1)=R0L(i)+1;
                RL0(i+1)=RL0(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00+Pk01_00 && P<Pk10_00+Pk01_00+Pk10_0L
                R0L(i+1)=R0L(i)-1;
                RLL(i+1)=RLL(i)+1;
                RL0(i+1)=RL0(i);
                R00(i+1)=R00(i);                
            elseif P>Pk10_00+Pk01_00+Pk10_0L && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0
                RL0(i+1)=RL0(i)-1;
                RLL(i+1)=RLL(i)+1;
                R0L(i+1)=R0L(i);
                R00(i+1)=R00(i);                
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0 && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00
                RL0(i+1)=RL0(i)-1;
                R00(i+1)=R00(i)+1;
                R0L(i+1)=R0L(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00 && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00
                R0L(i+1)=R0L(i)-1;
                R00(i+1)=R00(i)+1;
                RL0(i+1)=RL0(i);
                RLL(i+1)=RLL(i);
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00 && P<Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00+Pl10_LL
                RLL(i+1)=RLL(i)-1;
                R0L(i+1)=R0L(i)+1;
                RL0(i+1)=RL0(i);
                R00(i+1)=R00(i);
            elseif P>Pk10_00+Pk01_00+Pk10_0L+Pk01_L0+Pl10_00+Pl01_00+Pl10_LL
                RLL(i+1)=RLL(i)-1;
                RL0(i+1)=RL0(i)+1;
                R0L(i+1)=R0L(i);
                R00(i+1)=R00(i);
            end
            tita=(RL0(i)+R0L(i)+2*RLL(i))/(2*V)
            if tita_final-tita<0
                contador=contador+1;
            end
        end
        R00j(:,j)=interp1(t,R00,tj);R0Lj(:,j)=interp1(t,R0L,tj);RL0j(:,j)=interp1(t,RL0,tj);RLLj(:,j)=interp1(t,RLL,tj);
        clear R00 RL0 R0L RLL t
        n+m
        j
    end
    NC(m)=struct('RLL',RLLj,'tita',(2*RLLj+RL0j+R0Lj)/(2*V),'t',tj);
end
%% Guardo (LL la guardo como concentración en vez de Nº de partículas)
LL=LL/V;savefile = 'datastocastico.mat';save(savefile,'NC','IB','LL','V','-v7.3');
