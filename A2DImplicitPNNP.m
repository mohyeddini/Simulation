%AmirMohyeddini



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
clc;
clear;
close all;
phi=0.2;


leny=2000;
lenz=20;

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
qsc=zeros(Ny,Nx);
qsc(10,10)=-350;
qsc(15,20)=-360;

pl=6000;
pr=6000;

zarayeb=zeros(Nx*Ny,Nx*Ny);
sabet=zeros(Nx*Ny,1);
deltat=0.25;
time=45;
nt=time/deltat;
pressure=zeros(Ny,Nx,nt);
pressure(:,:,1)=6000;
bettac=1.127e-3;
alphac=5.615;
Ax=leny*lenz;

% dx=randi([100,400],[20,25]);
% dy=randi([100,400],[20,25]);
% load('matlabdx.mat');load('matlabdy.mat');



for n=2:nt
   n %#ok  
  for j=1:Ny
     for i=1:Nx
        index=i+(j-1)*Nx;
        vb=dx(j,i)*dy(j,i)*lenz;
   
        
        if i==1&&j==1
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i))))*((dx(j,i)+dx(j,i))/(dx(j,i)/perm(j,i)+dx(j,i)/perm(j,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
%             Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
            zarayeb(index,index+1)=Txr;
%             zarayeb(index,index-1)=Txl;
%             zarayeb(index,index-Nx)=Tyt;
            zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txr-Txl-Tyb-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb-Txl*pl;

            
        elseif i==Nx&&j==1
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i))))*((dx(j,i)+dx(j,i))/(dx(j,i)/perm(j,i)+dx(j,i)/perm(j,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
%             Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
%             zarayeb(index,index+1)=Txr;
            zarayeb(index,index-1)=Txl;
%             zarayeb(index,index-Nx)=Tyt;
            zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txl-Txr-Tyb-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb-Txr*pr;

            
        elseif j==Ny&&i==1     
            
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i))))*((dx(j,i)+dx(j,i))/(dx(j,i)/perm(j,i)+dx(j,i)/perm(j,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
%             Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
            zarayeb(index,index+1)=Txr;
%             zarayeb(index,index-1)=Txl;
            zarayeb(index,index-Nx)=Tyt;
%             zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txr-Txl-Tyt-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb-Txl*pl;


        elseif j==Ny&&i==Nx
            
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i))))*((dx(j,i)+dx(j,i))/(dx(j,i)/perm(j,i)+dx(j,i)/perm(j,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
%             Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
%             zarayeb(index,index+1)=Txr;
            zarayeb(index,index-1)=Txl;
            zarayeb(index,index-Nx)=Tyt;
%             zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txl-Txr-Tyt-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb-Txr*pr;
            
        
        elseif j==1
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
%             Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
            zarayeb(index,index+1)=Txr;
            zarayeb(index,index-1)=Txl;
%             zarayeb(index,index-Nx)=Tyt;
            zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txr-Txl-Tyb-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb;

            
        elseif j==Ny
            
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
%             Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
            zarayeb(index,index+1)=Txr;
            zarayeb(index,index-1)=Txl;
            zarayeb(index,index-Nx)=Tyt;
%             zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txr-Txl-Tyt-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb;

            

        elseif i==Nx
            
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i))))*((dx(j,i)+dx(j,i))/(dx(j,i)/perm(j,i)+dx(j,i)/perm(j,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
%             zarayeb(index,index+1)=Txr;
            zarayeb(index,index-1)=Txl;
            zarayeb(index,index-Nx)=Tyt;
            zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txl-Txr-Tyt-Tyb-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb-Txr*pr;


        elseif i==1    
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i))))*((dx(j,i)+dx(j,i))/(dx(j,i)/perm(j,i)+dx(j,i)/perm(j,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
            zarayeb(index,index+1)=Txr;
%             zarayeb(index,index-1)=Txl;
            zarayeb(index,index-Nx)=Tyt;
            zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txr-Txl-Tyt-Tyb-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb-Txl*pl;

        else
            
            Txr=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i+1))))*((dx(j,i)+dx(j,i+1))/(dx(j,i)/perm(j,i)+dx(j,i+1)/perm(j,i+1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i+1,n-1)))));
            Txl=bettac*(1/dx(j,i)/(1/2*(dx(j,i)+dx(j,i-1))))*((dx(j,i)+dx(j,i-1))/(dx(j,i)/perm(j,i)+dx(j,i-1)/perm(j,i-1)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i-1,n-1)))));
            Tyt=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j-1,i))))*((dy(j,i)+dy(j-1,i))/(dy(j,i)/perm(j,i)+dy(j-1,i)/perm(j-1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j-1,i,n-1)))));
            Tyb=bettac*(1/dy(j,i)/(1/2*(dy(j,i)+dy(j+1,i))))*((dy(j,i)+dy(j+1,i))/(dy(j,i)/perm(j,i)+dy(j+1,i)/perm(j+1,i)))*(1/(mu(1/2*(pressure(j,i,n-1)+pressure(j+1,i,n-1)))*Bo(1/2*(pressure(j,i,n-1)+pressure(j,i,n-1)))));
           
            
            zarayeb(index,index+1)=Txr;
            zarayeb(index,index-1)=Txl;
            zarayeb(index,index-Nx)=Tyt;
            zarayeb(index,index+Nx)=Tyb;
            zarayeb(index,index)=-Txr-Txl-Tyt-Tyb-phi*co/(alphac*deltat);
            sabet(index,1)=-phi*co/(alphac*deltat)*pressure(j,i,n-1)-qsc(j,i)/vb;
        end
    end
     
  end
   
  
  temp=linsolve(zarayeb,sabet);
  for j=1:Ny
      for i=1:Nx
      index=i+(j-1)*Nx;
      pressure(j,i,n)=temp(index);
      end
  end

end
% for n=2:nt
% pcolor(pressure(:,:,n))
% drawnow
% pause(0.025)
% end

for n=2:nt
surf(pressure(:,:,n))
 zlim([4000,6000])
drawnow
pause(0.025)
end
