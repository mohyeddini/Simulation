%Amir Mohyeddini
%EXPLICIT
%noflow______noflow
%else is important




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
clc;
clear;
close all;
N=5;
phi=0.18;

deltax=zeros(1,N);
for i=1:N
   deltax(1,i)=1000; 
end
dy=1000;
dz=75;
% permx=15;
permx=zeros(1,N);
for i=1:N
   permx(1,i)=15; 
    
end
% Bo=1;
co=3.5e-6;
Bo=@(p) 1;%1/(1+co*(p-4200));
muob=3.5;
pb=2635;
mu=@(p) 10;%muob*(p/pb)^(2.6*(p^1.187)*exp(-11.513+(-8.98e-5*p)));


deltat=15;
time=360;
nt=time/deltat;
pressure=zeros(nt,N);
pint=6000;
pressure(1,:)=pint;
betac=1.127e-3;
alphac=5.615;
Ax=dy*dz;
vb=Ax.*deltax;
q=zeros(1,N);
q(4)=-150;
pr=6000;
pl=6000;

for n=2:nt

   for i=1:N 

        if i==1
         
         right=(betac*Ax*(deltax(1,i)+deltax(1,i+1))/(deltax(1,i)/permx(1,i)+deltax(1,i+1)/permx(1,i+1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*(1/2*(deltax(1,i)+deltax(1,i+1))));
         left=(betac*Ax*(deltax(1,i)+deltax(1,i))/(deltax(1,i)/permx(1,i)+deltax(1,i)/permx(1,i)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i)))*(1/2*(deltax(1,i)+deltax(1,i))));

         pressure(n,i)=pressure(n-1,i)+(alphac*Bo(pressure(n-1,i))*deltat)/(vb(1,i)*phi*co)*q(1,i)+(alphac*Bo(pressure(n-1,i))*deltat)/(vb(1,i)*phi*co)*((right*pressure(n-1,i+1)-(right+left)*(pressure(n-1,i)))+(left*pl));   

            
        elseif i==N
          
         right=(betac*Ax*(deltax(1,i)+deltax(1,i))/(deltax(1,i)/permx(1,i)+deltax(1,i)/permx(1,i)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i)))*(1/2*(deltax(1,i)+deltax(1,i))));
         left=(betac*Ax*(deltax(1,i)+deltax(1,i-1))/(deltax(1,i)/permx(1,i)+deltax(1,i-1)/permx(1,i-1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*(1/2*(deltax(1,i)+deltax(1,i-1))));

        pressure(n,i)=pressure(n-1,i)+(alphac*Bo(pressure(n-1,i))*deltat)/(vb(1,i)*phi*co)*q(1,i)+(alphac*Bo(pressure(n-1,i))*deltat)/(vb(1,i)*phi*co)*((right*pr-(right+left)*(pressure(n-1,i)))+(left*pressure(n-1,i-1)));   
  
            
        else
         right=(betac*Ax*(deltax(1,i)+deltax(1,i+1))/(deltax(1,i)/permx(1,i)+deltax(1,i+1)/permx(1,i+1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*(1/2*(deltax(1,i)+deltax(1,i+1))));
         left=(betac*Ax*(deltax(1,i)+deltax(1,i-1))/(deltax(1,i)/permx(1,i)+deltax(1,i-1)/permx(1,i-1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*(1/2*(deltax(1,i)+deltax(1,i-1))));

        pressure(n,i)=pressure(n-1,i)+(alphac*Bo(pressure(n-1,i))*deltat)/(vb(1,i)*phi*co)*q(1,i)+(alphac*Bo(pressure(n-1,i))*deltat)/(vb(1,i)*phi*co)*((right*pressure(n-1,i+1)-(right+left)*(pressure(n-1,i)))+(left*pressure(n-1,i-1)));   


        end
    
    end
   
end

for n=1:nt
   plot(pressure(n,:))
   hold on
end


