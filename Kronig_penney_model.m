clc;
clear all;
close all;
%Physical Constants
hbar = 1.055e-34; %hbar in units of J.s
m0 = 9.109e-31; %mass of an electron in kg
eV = 1.602e-19; %Electron Volt in J
%Problem Specific Constants
a = 0.1e-9; %a,b from Kronig-Penny Model in m
b = 0.5e-9;
d=a+b;
U0 = 3*eV; %Height of periodic potenial square wave in J.
dE = 0.00001*eV; %Set dE to a small increment
E_value = (0:dE:20*eV);
E_value_eV = (0:dE/eV:20);

zita = zeros(size(E_value));
for i=1:length(zita)
    zita(i)=E_value(i)/U0;
end
F=zeros(size(E_value));
alpha_0=sqrt((2*m0*U0)/hbar^2);

for i=1:length(zita)
    if zita(i)<1
        F(i)=((1-2*zita(i))/(2*sqrt(zita(i)*(1-zita(i)))))*sin(alpha_0*a*sqrt(zita(i)))*sinh(alpha_0*b*sqrt(1-zita(i)))+cos(alpha_0*a*sqrt(zita(i)))*cosh(alpha_0*b*sqrt(1-zita(i)));
    else
        F(i)=((1-2*zita(i))/(2*sqrt(zita(i)*(zita(i)-1))))*sin(alpha_0*a*sqrt(zita(i)))*sin(alpha_0*b*sqrt(zita(i)-1))+cos(alpha_0*a*sqrt(zita(i)))*cos(alpha_0*b*sqrt(zita(i)-1));
    end
    k(i)=acos(F(i))/d;
end
E_1=zeros(1,10);
One_ve=zeros(size(zita));
One_minus_ve=zeros(size(zita));
for i=1:length(zita)
    One_ve(i)=1;
    One_minus_ve(i)=-1;
end
tolerance=0.00001;
loop_count=1;
for i=1:length(zita)
    if abs(F(i)-One_minus_ve(i))<tolerance | abs(F(i)-One_ve(i))<tolerance
        E_1(loop_count)=zita(i);
        loop_count=loop_count+1;
    end
end

figure(1)
plot(k, E_value_eV);
xlabel('k','fontweight','bold','fontsize',16),ylabel('E (eV)','fontweight','bold','fontsize',16)

figure(2)
plot(zita, F);
hold on;

plot(zita,One_ve);
hold on;

plot(zita,One_minus_ve);
xlabel('E/U0','fontweight','bold','fontsize',16),ylabel('left side function','fontweight','bold','fontsize',16)
