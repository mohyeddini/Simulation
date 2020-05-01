%AmirMohyeddini




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
clc;
clear;
close all;
phi=0.2;

lenx=2000;
leny=2000;
lenz=20;
pbh=3600;
load('matlab.mat');

co=3.5e-6;
Bo=@(p) 1;%1/(1+co*(p-4200));
muob=3.5;
pb=2635;
mu=@(p) muob*(p/pb)^(2.6*(p^1.187)*exp(-11.513+(-8.98e-5*p)));


Nx=25;
Ny=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dx=lenx/Nx;

dx=randi([100,300],[Ny,Nx]);
dy=randi([100,300],[Ny,Nx]);



% qsc=zeros(Ny,Nx);
% qsc(5,5)=-50;
% qsc(15,20)=-75;

deltat=0.01;
time=100;
nt=time/deltat;
pressure=zeros(Ny,Nx,nt);
pressure(:,:,1)=6000;
bettac=1.127e-3;
alphac=5.615;
Ax=leny*lenz;


rw=0.25;
for n=2:nt
 
    n %#ok
  for i=1:Nx
     for j=1:Ny
        re=sqrt(dx(j,i)*dy(j,i)/pi);
        wc=-bettac*2*pi*perm(j,i)*lenz/(mu(pressure(j,i,n-1))*Bo(pressure(j,i,n-1))*log(re/rw));
        vb=dx(j,i)*dy(j,i)*lenz;
        if i==1&&j==1
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i+1)+dx(j,i))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
%             Txl=bettac*(1/dx^2)*(2/(1/perm(j,i)+1/perm(j,i-1)))*(1/(mu*Bo(pressure(j,i,n-1))));
%             Tyt=bettac*(1/dy^2)*(2/(1/perm(j,i)+1/perm(j+1,i)))*(1/(mu*Bo(pressure(j,i,n-1))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=Txr*(pressure(j,i+1,n-1)-pressure(j,i,n-1));
            sharp=-Tyb*(pressure(j,i,n-1)-pressure(j+1,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));  
            
        elseif i==Nx&&j==1
%             Txr=bettac*(1/dx^2)*(2/(1/perm(j,i)+1/perm(j,i+1)))*(1/(mu*Bo(pressure(j,i,n-1))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
%             Tyt=bettac*(1/dy^2)*(2/(1/perm(j,i)+1/perm(j+1,i)))*(1/(mu*Bo(pressure(j,i,n-1))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=-Txl*(pressure(j,i,n-1)-pressure(j,i-1,n-1));
            sharp=-Tyb*(pressure(j,i,n-1)-pressure(j+1,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        elseif j==Ny&&i==1     
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
%             Txl=bettac*(1/(1/2*(dx(j,i)+dx(j,i-1)))^2)*(dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1))*(1/(mu*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
%             Tyb=bettac*(1/(1/2*(dy(j,i)+dy(j+1,i)))^2)*((dy(j,i)+dy(j,i+1))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=Txr*(pressure(j,i+1,n-1)-pressure(j,i,n-1));
            sharp=Tyt*(pressure(j-1,i,n-1)-pressure(j,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        elseif j==Ny&&i==Nx
%             Txr=bettac*(1/dx^2)*(2/(1/perm(j,i)+1/perm(j,i+1)))*(1/(mu*Bo(pressure(j,i,n-1))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
%             Tyb=bettac*(1/(1/2*(dy(j,i)+dy(j+1,i)))^2)*((dy(j,i)+dy(j,i+1))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=-Txl*(pressure(j,i,n-1)-pressure(j,i-1,n-1));
            sharp=Tyt*(pressure(j-1,i,n-1)-pressure(j,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        elseif j==1
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
%             Tyt=bettac*(1/(1/2*(dy(j,i)+dy(j-1,i))^2)*(dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=Txr*(pressure(j,i+1,n-1)-pressure(j,i,n-1))-Txl*(pressure(j,i,n-1)-pressure(j,i-1,n-1));
            sharp=-Tyb*(pressure(j,i,n-1)-pressure(j+1,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        elseif j==Ny
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
%             Tyb=bettac*(1/(1/2*(dy(j,i)+dy(j+1,i)))^2)*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=Txr*(pressure(j,i+1,n-1)-pressure(j,i,n-1))-Txl*(pressure(j,i,n-1)-pressure(j,i-1,n-1));
            sharp=Tyt*(pressure(j-1,i,n-1)-pressure(j,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        elseif i==Nx
%             Txr=bettac*(1/(1/2*(dx(j,i)+dx(j,i+1)))^2)*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=-Txl*(pressure(j,i,n-1)-pressure(j,i-1,n-1));
            sharp=Tyt*(pressure(j-1,i,n-1)-pressure(j,i,n-1))-Tyb*(pressure(j,i,n-1)-pressure(j+1,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        elseif i==1    
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
%             Txl=bettac*(1/(1/2*(dx(j,i)+dx(j,i-1)))^2)*(dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1))*(1/(mu*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=Txr*(pressure(j,i+1,n-1)-pressure(j,i,n-1));
            sharp=Tyt*(pressure(j-1,i,n-1)-pressure(j,i,n-1))-Tyb*(pressure(j,i,n-1)-pressure(j+1,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        else
        
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))));
            star=Txr*(pressure(j,i+1,n-1)-pressure(j,i,n-1))-Txl*(pressure(j,i,n-1)-pressure(j,i-1,n-1));
            sharp=Tyt*(pressure(j-1,i,n-1)-pressure(j,i,n-1))-Tyb*(pressure(j,i,n-1)-pressure(j+1,i,n-1));
            pressure(j,i,n)=pressure(j,i,n-1)+star*(alphac*1*deltat)/(phi*co)+sharp*(alphac*Bo(pressure(j,i,n-1))*deltat/(phi*co));
        
        if i==10&&j==10
            pressure(j,i,n)=pressure(j,i,n)+wc*(pressure(j,i,n-1)-pbh)*(alphac*Bo(pressure(j,i,n-1))*deltat)/(vb*phi*co);
        
        end
        
        end
    end
     
  end

end

% for n=2:nt
% pcolor(pressure(:,:,n))
% drawnow
% pause(0.025)
% end

%  plot section

for n=2:nt
surf(pressure(:,:,n))
%  zlim([4000,6000])
drawnow
% pause(0.025)
end
