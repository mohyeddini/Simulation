%AmirMohyeddini

%comment
%press|_______|noflow
                                %point distribute
%zarib
%sabet
%zarib=zeros(N-1,N-1);
%sabet=zeros(N-1,1);
%pressure(:,1)=pl;
%if i==1 sabet -left*pl
%press_________noflow

%                                               point distribution
%pressure(n-1,i+1)&&-qsc(i+1) important
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

N=5;
phi=0.18;

deltax=zeros(1,N);
for i=1:N
   deltax(1,i)=1000; 
end
deltay=1000;
deltaz=75;
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


% deltax=lenx/N;
zarib=zeros(N-2,N-2);
sabet=zeros(N-2,1);
deltat=15;
time=365;
pl=6000;
pr=6000;
nt=floor(time/deltat);
pressure=zeros(nt,N);
pressure(1,:)=6000;
pressure(:,1)=pl;
pressure(:,end)=6000;
betac=1.127e-3;
alphac=5.615;
Ax=deltay*deltaz;
vb=Ax.*deltax;

qsc=zeros(1,N);
qsc(4)=-150;

% qsc(67)=-10;
for n=2:nt
   for i=1:N-2
        if i==1
                     
         right=(betac*Ax*(deltax(1,i+1)+deltax(1,i+2))/(deltax(1,i+1)/permx(1,i+1)+deltax(1,i+2)/permx(1,i+2)))/(mu(1/2*(pressure(n-1,i+1)+pressure(n-1,i+2)))*Bo(1/2*(pressure(n-1,i+1)+pressure(n-1,i+2)))*(1/2*(deltax(1,i+1)+deltax(1,i+2))));
         left=(betac*Ax*(deltax(1,i+1)+deltax(1,i))/(deltax(1,i+1)/permx(1,i+1)+deltax(1,i)/permx(1,i)))/(mu(1/2*(pressure(n-1,i+1)+pressure(n-1,i)))*Bo(1/2*(pressure(n-1,i+1)+pressure(n-1,i)))*(1/2*(deltax(1,i+1)+deltax(1,i))));
         self=-2*(betac*Ax*permx(1,i+1))/(mu(pressure(n-1,i+1))*Bo(pressure(n-1,i+1))*(deltax(1,i+1)))-(vb(1,i+1)*phi*co)/(alphac*Bo(pressure(n-1,i+1))*deltat);  
           zarib(i,i+1)=right;
           zarib(i,i)=self;
%            zarib(i,i-1)=left;
           sabet(i,1)=-qsc(i+1)-(vb(1,i+1)*phi*co)/(alphac*Bo(pressure(n-1,i+1))*deltat)*pressure(n-1,i+1)-pl*left;

      elseif i==N-2
                   
         right=(betac*Ax*(deltax(1,i+1)+deltax(1,i+2))/(deltax(1,i+1)/permx(1,i+1)+deltax(1,i+2)/permx(1,i+2)))/(mu(1/2*(pressure(n-1,i+1)+pressure(n-1,i+2)))*Bo(1/2*(pressure(n-1,i+1)+pressure(n-1,i+2)))*(1/2*(deltax(1,i+1)+deltax(1,i+2))));
         left=(betac*Ax*(deltax(1,i+1)+deltax(1,i))/(deltax(1,i+1)/permx(1,i+1)+deltax(1,i)/permx(1,i)))/(mu(1/2*(pressure(n-1,i+1)+pressure(n-1,i)))*Bo(1/2*(pressure(n-1,i+1)+pressure(n-1,i)))*(1/2*(deltax(1,i+1)+deltax(1,i))));
         self=-2*(betac*Ax*permx(1,i+1))/(mu(pressure(n-1,i+1))*Bo(pressure(n-1,i+1))*(deltax(1,i+1)))-(vb(1,i+1)*phi*co)/(alphac*Bo(pressure(n-1,i+1))*deltat);  
%             zarib(i,i+1)=right;
           zarib(i,i)=self;
           zarib(i,i-1)=left;
           sabet(i,1)=-qsc(i+1)-(vb(1,i+1)*phi*co)/(alphac*Bo(pressure(n-1,i+1))*deltat)*pressure(n-1,i+1)-pr*right;
        else
            
         right=(betac*Ax*(deltax(1,i+1)+deltax(1,i+2))/(deltax(1,i+1)/permx(1,i+1)+deltax(1,i+2)/permx(1,i+2)))/(mu(1/2*(pressure(n-1,i+1)+pressure(n-1,i+2)))*Bo(1/2*(pressure(n-1,i+1)+pressure(n-1,i+2)))*(1/2*(deltax(1,i+1)+deltax(1,i+2))));
         left=(betac*Ax*(deltax(1,i+1)+deltax(1,i))/(deltax(1,i+1)/permx(1,i+1)+deltax(1,i)/permx(1,i)))/(mu(1/2*(pressure(n-1,i+1)+pressure(n-1,i)))*Bo(1/2*(pressure(n-1,i+1)+pressure(n-1,i)))*(1/2*(deltax(1,i+1)+deltax(1,i))));
         self=-2*(betac*Ax*permx(1,i+1))/(mu(pressure(n-1,i+1))*Bo(pressure(n-1,i+1))*(deltax(1,i+1)))-(vb(1,i+1)*phi*co)/(alphac*Bo(pressure(n-1,i+1))*deltat);  
           zarib(i,i+1)=right;
           zarib(i,i)=self;
           zarib(i,i-1)=left;
           sabet(i,1)=-qsc(i+1)-(vb(1,i+1)*phi*co)/(alphac*Bo(pressure(n-1,i+1))*deltat)*pressure(n-1,i+1);
           
       end
       
      
   end
   %x=zarayeb^-1*sabet;
x=linsolve(zarib,sabet);
pressure(n,2:end-1)=x';
% pressure(:,1)=6000;
    
    
end
% for n=2:nt
%     (pressure(n,5))-(pressure(n-1,5));  %#ok
% end

for n=1:nt
   plot(pressure(n,:))
   hold on
end


% figure
% plot(pressure(:,50));

