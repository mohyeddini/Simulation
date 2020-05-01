%Amir Mohyeddini

%comment
    
%noflow|_____|noflow

    %deltax
    %zarib
    %if i==1 left*pi
    %if i==N right*pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

zarib=zeros(N,N);%zarib
sabet=zeros(N,1);%sabet
deltat=1;
time=45;
nt=time/deltat;
pressure=zeros(nt,N);
pressure(1,:)=6000;
betac=1.127e-3;
alphac=5.615;
Ax=deltay*deltaz;
vb=Ax.*deltax;
pl=6000;
pr=6000;

qsc=zeros(1,N);%q
qsc(4)=-5;
qsc(2)=-10;

for n=2:nt
   for i=1:N
       if i==1
           
        right=(betac*Ax*(deltax(1,i)+deltax(1,i+1))/(deltax(1,i)/permx(1,i)+deltax(1,i+1)/permx(1,i+1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*(1/2*(deltax(1,i)+deltax(1,i+1))));
         left=(betac*Ax*(deltax(1,i)+deltax(1,i))/(deltax(1,i)/permx(1,i)+deltax(1,i)/permx(1,i)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i)))*(1/2*(deltax(1,i)+deltax(1,i))));
         self=-2*(betac*Ax*permx(1,i))/(mu(pressure(n-1,i))*Bo(pressure(n-1,i))*(deltax(1,i)))-(vb(1,i)*phi*co)/(alphac*Bo(pressure(n-1,i))*deltat); 
         zarib(i,i+1)=right;
           zarib(i,i)=self;
%            zarib(i,i-1)=left;
           sabet(i,1)=-qsc(i)-(vb(1,i)*phi*co)/(alphac*Bo(pressure(n-1,i))*deltat)*pressure(n-1,i)-pl*left;           
           
       elseif i==N
        right=(betac*Ax*(deltax(1,i)+deltax(1,i))/(deltax(1,i)/permx(1,i)+deltax(1,i)/permx(1,i)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i)))*(1/2*(deltax(1,i)+deltax(1,i))));
         left=(betac*Ax*(deltax(1,i)+deltax(1,i-1))/(deltax(1,i)/permx(1,i)+deltax(1,i-1)/permx(1,i-1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*(1/2*(deltax(1,i)+deltax(1,i-1))));
         self=-2*(betac*Ax*permx(1,i))/(mu(pressure(n-1,i))*Bo(pressure(n-1,i))*(deltax(1,i)))-(vb(1,i)*phi*co)/(alphac*Bo(pressure(n-1,i))*deltat); 
%          zarib(i,i+1)=right;
           zarib(i,i)=self;
           zarib(i,i-1)=left;
           sabet(i,1)=-qsc(i)-(vb(1,i)*phi*co)/(alphac*Bo(pressure(n-1,i))*deltat)*pressure(n-1,i)-pr*right;           
       else
         right=(betac*Ax*(deltax(1,i)+deltax(1,i+1))/(deltax(1,i)/permx(1,i)+deltax(1,i+1)/permx(1,i+1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i+1)))*(1/2*(deltax(1,i)+deltax(1,i+1))));
         left=(betac*Ax*(deltax(1,i)+deltax(1,i-1))/(deltax(1,i)/permx(1,i)+deltax(1,i-1)/permx(1,i-1)))/(mu(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*Bo(1/2*(pressure(n-1,i)+pressure(n-1,i-1)))*(1/2*(deltax(1,i)+deltax(1,i-1))));
         self=-2*(betac*Ax*permx(1,i))/(mu(pressure(n-1,i))*Bo(pressure(n-1,i))*(deltax(1,i)))-(vb(1,i)*phi*co)/(alphac*Bo(pressure(n-1,i))*deltat); 
         zarib(i,i+1)=right;
           zarib(i,i)=self;
           zarib(i,i-1)=left;
           sabet(i,1)=-qsc(i)-(vb(1,i)*phi*co)/(alphac*Bo(pressure(n-1,i))*deltat)*pressure(n-1,i);           
       end
       
      
   end
   %x=zarayeb^-1*sabet;
x=linsolve(zarib,sabet);
pressure(n,:)=x';

    
    
end

% for n=2:nt
%     (pressure(n,50))-(pressure(n-1,50));  %#ok
% end

for n=1:nt
   plot(pressure(n,:))
   hold on
end


% figure
% plot(pressure(:,50));

