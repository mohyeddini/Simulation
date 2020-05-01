%AmirMohyeddini




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1

clc
clear
close all;

p1 =-4.167;
p2=12.5;
p3=-15.46;
p4=7.125;

pc_f=@(s) (p1*s.^3 + p2*s.^2 + p3*s+p4)*6894.7573;%pa
dpc_dsw=@(s) (3.*p1*s.^2 + 2.*p2*s + p3)*6894.7573;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lenx=200*0.3048;
leny=200*0.3048;
lenz=50*0.3048;%m
dx1=10*0.3048;
dy1=10*0.3048;
dz=lenz;%m
time_total=1*100*24*60*60;%sec
dt=1*24*60*60;%sec
Nt=round(time_total/dt);
Nx=lenx/dx1;
Ny=leny/dy1;


dx=zeros(1,Nx);
dy=zeros(1,Nx);
dx(1,:)=dx1;
dy(1,:)=dy1;
%  Nx=20;
%  Ny=20;

% dx1=randi([10,20],[Ny,Nx]);
% dy1=randi([10,20],[Ny,Nx]);
%  dx=dx1.*0.3048;
%  dy=dy1.*0.3048;

rw=0.3*0.3048;%meter
re=sqrt((dx.*dy)./pi);%m

perm=xlsread('perm.xlsx');
perm=perm(1,1:Nx)*0.001*9.869233e-13;%m2
% perm=ones(Ny,Nx)*0.5*9.869233e-13;

bo0=1;
bw0=1;

Cw=1/6894.7573*1e-5;
Co=1/6894.7573*1e-5;
Cr=1/6894.7573*1e-8;

fi=xlsread('fi.xlsx');
fi=fi(1,1:Nx);

Bo=@(p) 1; %1./(1+Co*(p-3200*6894.7573));
Bw=@(p) 1 ;%1./(1+Cw*(p-3200*6894.7573));

mu_o=@(p) 0.005;%
mu_w=@(p) 0.001;%           %pa.s

% muob=3.5;
% pbo=2635;
% mu_o=@(p) (muob*(p/pbo)^(2.6*(p^1.187)*exp(-11.513+(-8.98e-5*p))))/10^3;
% 
% muwb=3.5;
% pbw=2635;
% mu_w=@(p) (muob*(p/pbw)^(2.6*(p^1.187)*exp(-11.513+(-8.98e-5*p))))/10^3;

% initial conditions
swi=0.2;
swc=0.2;
sor=0.2;
init_p=5000*6894.7573;%pa
pressure=zeros(Nt,Nx);
sat=zeros(Nt,Nx);
pressure(1,:)=init_p;
sat(1,:)=swi;
zarayeb=zeros(Nx,Nx);
sabet=zeros(Nx,1);
trans=@(kr1,kr2,mu1,mu2,B1,B2,pre1,pre2,k1,k2,d1,d2) 2*((d1*kr1/mu1/B1+d2*kr2/mu2/B2)/(d1+d2))/(d1)/(d2/k2+d1/k1);
q_inj=5*5.615*(0.3048)^3/86400;
pbh=3600*6894.7573;



for n=2:Nt
    n %#ok
    
        for i=1:Nx
            
            landa_o=rel_perm(sat(n-1,i),2)/mu_o(pressure(n-1,i))/Bo(pressure(n-1,i));
            landa_w=rel_perm(sat(n-1,i),1)/mu_w(pressure(n-1,i))/Bw(pressure(n-1,i));
            cpooi=fi(1,i)*(1-sat(n-1,i))/dt*(Cr/Bo(pressure(n-1,i))+Co);
            cswoi=-fi(1,i)/dt/Bo(pressure(n-1,i));
            cpowi=fi(1,i)*(sat(n-1,i))/dt*(Cr/Bw(pressure(n-1,i))+Cw);
            cswwi=fi(1,i)/dt/Bw(pressure(n-1,i))-dpc_dsw(sat(n-1,i))*cpowi;
            alpha=-cswwi/cswoi;
            
            
              if i==1
               trans_x_o_r=trans(rel_perm(sat(n-1,i),2),rel_perm(sat(n-1,i+1),2),mu_o(pressure(n-1,i)),mu_o(pressure(n-1,i+1)),Bo(pressure(n-1,i)),Bo(pressure(n-1,i+1)),pressure(n-1,i),pressure(n-1,i+1),perm(1,i),perm(1,i+1),dx(1,i),dx(1,i+1));
               trans_x_w_r=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i+1),1),mu_w(pressure(n-1,i)),mu_w(pressure(n-1,i+1)),Bw(pressure(n-1,i)),Bw(pressure(n-1,i+1)),pressure(n-1,i),pressure(n-1,i+1),perm(1,i),perm(1,i+1),dx(1,i),dx(1,i+1));
              A=-(trans_x_o_r+cpooi)-alpha*(trans_x_w_r+cpowi);
              B=trans_x_o_r+alpha*trans_x_w_r;
              F=-(cpooi+alpha*cpowi)*pressure(n-1,i)+alpha*(trans_x_w_r*(pc_f(sat(n-1,i+1))-pc_f(sat(n-1,i))));
                F=F-alpha*abs(q_inj)/dx(1,i)/dy(1,i)/dz;
                zarayeb(i,i+1)=B;
                zarayeb(i,i)=A;
                sabet(i,1)=F;
            elseif i==Nx
                 trans_x_o_l=trans(rel_perm(sat(n-1,i),2),rel_perm(sat(n-1,i-1),2),mu_o(pressure(n-1,i)),mu_o(pressure(n-1,i-1)),Bo(pressure(n-1,i)),Bo(pressure(n-1,i-1)),pressure(n-1,i),pressure(n-1,i-1),perm(1,i),perm(1,i-1),dx(1,i),dx(1,i-1));
                 trans_x_w_l=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i-1),1),mu_w(pressure(n-1,i)),mu_w(pressure(n-1,i-1)),Bw(pressure(n-1,i)),Bw(pressure(n-1,i-1)),pressure(n-1,i),pressure(n-1,i-1),perm(1,i),perm(1,i-1),dx(1,i),dx(1,i-1));
                A=-(trans_x_o_l+cpooi)-alpha*(trans_x_w_l+cpowi);
                C=trans_x_o_l+alpha*trans_x_w_l;
                 F=-(cpooi+alpha*cpowi)*pressure(n-1,i)+alpha*(trans_x_w_l*(pc_f(sat(n-1,i-1))-pc_f(sat(n-1,i))));
                 WC=2*pi*perm(1,i)*dz/log(re(1,i)/rw);            
                 A=A-WC/dx(1,i)/dy(1,i)/dz*landa_o-alpha*WC/dx(1,i)/dy(1,i)/dz*landa_w;
                F=F-WC/dx(1,i)/dy(1,i)/dz*landa_o*pbh-alpha*WC/dx(1,i)/dy(1,i)/dz*landa_w*pbh;
                zarayeb(i,i-1)=C;
                zarayeb(i,i)=A;
                sabet(i,1)=F;
              else
                 trans_x_o_r=trans(rel_perm(sat(n-1,i),2),rel_perm(sat(n-1,i+1),2),mu_o(pressure(n-1,i)),mu_o(pressure(n-1,i+1)),Bo(pressure(n-1,i)),Bo(pressure(n-1,i+1)),pressure(n-1,i),pressure(n-1,i+1),perm(1,i),perm(1,i+1),dx(1,i),dx(1,i+1));
                 trans_x_w_r=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i+1),1),mu_w(pressure(n-1,i)),mu_w(pressure(n-1,i+1)),Bw(pressure(n-1,i)),Bw(pressure(n-1,i+1)),pressure(n-1,i),pressure(n-1,i+1),perm(1,i),perm(1,i+1),dx(1,i),dx(1,i+1));
                 trans_x_o_l=trans(rel_perm(sat(n-1,i),2),rel_perm(sat(n-1,i-1),2),mu_o(pressure(n-1,i)),mu_o(pressure(n-1,i-1)),Bo(pressure(n-1,i)),Bo(pressure(n-1,i-1)),pressure(n-1,i),pressure(n-1,i-1),perm(1,i),perm(1,i-1),dx(1,i),dx(1,i-1));
                 trans_x_w_l=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i-1),1),mu_w(pressure(n-1,i)),mu_w(pressure(n-1,i-1)),Bw(pressure(n-1,i)),Bw(pressure(n-1,i-1)),pressure(n-1,i),pressure(n-1,i-1),perm(1,i),perm(1,i-1),dx(1,i),dx(1,i-1));
                A=-(trans_x_o_r+trans_x_o_l+cpooi)-alpha*(trans_x_w_r+trans_x_w_l+cpowi);
                B=trans_x_o_r+alpha*trans_x_w_r;
                C=trans_x_o_l+alpha*trans_x_w_l;
                F=-(cpooi+alpha*cpowi)*pressure(n-1,i)+alpha*(trans_x_w_r*(pc_f(sat(n-1,i+1))-pc_f(sat(n-1,i)))+trans_x_w_l*(pc_f(sat(n-1,i-1))-pc_f(sat(n-1,i))));
                
                zarayeb(i,i+1)=B;
                zarayeb(i,i-1)=C;
                zarayeb(i,i)=A;
                sabet(i,1)=F;
              end
              
              
        
        end
        
    javab=linsolve(zarayeb,sabet);
    
        
            
            pressure(n,:)=javab';
        
   
    
        for i=1:Nx
             if i==1
              trans_x_w_r=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i+1),1),mu_w(pressure(n,i)),mu_w(pressure(n,i+1)),Bw(pressure(n,i)),Bw(pressure(n,i+1)),pressure(n,i),pressure(n,i+1),perm(1,i),perm(1,i+1),dx(1,i),dx(1,i+1));
               sat(n,i)=sat(n-1,i)+1/cswwi*(trans_x_w_r*(pressure(n,i+1)-pc_f(sat(n-1,i+1))-(pressure(n,i)-pc_f(sat(n-1,i))))-cpowi*(pressure(n,i)-pressure(n-1,i)));
               sat(n,i)=sat(n,i)+1/cswwi/dx(1,i)/dy(1,i)/dz*q_inj; 
            elseif i==Nx
              trans_x_w_l=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i-1),1),mu_w(pressure(n,i)),mu_w(pressure(n,i-1)),Bw(pressure(n,i)),Bw(pressure(n,i-1)),pressure(n,i),pressure(n,i-1),perm(1,i),perm(1,i-1),dx(1,i),dx(1,i-1));
               sat(n,i)=sat(n-1,i)+1/cswwi*(trans_x_w_l*(pressure(n,i-1)-pc_f(sat(n-1,i-1))-(pressure(n,i)-pc_f(sat(n-1,i))))-cpowi*(pressure(n,i)-pressure(n-1,i)));
                
                 WC=2*pi*perm(1,i)*dz/log(re(1,i)/rw);            
                 sat(n,i)=sat(n,i)-1/cswwi*WC/dx(1,i)/dy(1,i)/dz*landa_w*(pressure(n,i)-pbh);
     
            else
              trans_x_w_r=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i+1),1),mu_w(pressure(n,i)),mu_w(pressure(n,i+1)),Bw(pressure(n,i)),Bw(pressure(n,i+1)),pressure(n,i),pressure(n,i+1),perm(1,i),perm(1,i+1),dx(1,i),dx(1,i+1));
              trans_x_w_l=trans(rel_perm(sat(n-1,i),1),rel_perm(sat(n-1,i-1),1),mu_w(pressure(n,i)),mu_w(pressure(n,i-1)),Bw(pressure(n,i)),Bw(pressure(n,i-1)),pressure(n,i),pressure(n,i-1),perm(1,i),perm(1,i-1),dx(1,i),dx(1,i-1));
             
                sat(n,i)=sat(n-1,i)+1/cswwi*(trans_x_w_r*(pressure(n,i+1)-pc_f(sat(n-1,i+1))-(pressure(n,i)-pc_f(sat(n-1,i))))+trans_x_w_l*(pressure(n,i-1)-pc_f(sat(n-1,i))-(pressure(n,i)-pc_f(sat(n-1,i))))-cpowi*(pressure(n,i)-pressure(n-1,i)));
            end
        end
    
end


for n=1:Nt
   plot(pressure(n,:))
   hold on
end

figure
for n=1:Nt
   plot(sat(n,:))
   hold on
end