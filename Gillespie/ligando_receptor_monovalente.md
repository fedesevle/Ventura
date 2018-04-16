% Modelo Ligando-Receptor monovalente estocástico que usa el algoritmo de Gillespie
% Lenguaje: MATLAB
clear all;close all
n=1000;
k=1;l=1;R0=100;RL=0;L=1000;t=0;V=100;
for i=1:n
    dR(i)=(k*R0(i)*L/V+l*RL(i)); %Probabilidad de que haya una reacción
    dt=exprnd(1/dR(i));
    t(i+1)=t(i)+dt;
    Pk=k*R0(i)*L/V/dR(i);
    Pl=1-Pk;
    P=rand();
    if P<Pk
        R0(i+1)=R0(i)-1;
        RL(i+1)=RL(i)+1;
    else
        RL(i+1)=RL(i)-1;
        R0(i+1)=R0(i)+1;
    end
end

plot(t,RL)
