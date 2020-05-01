%AmirMohyeddini






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
clc;
clear;
close all;

lenx=200*0.3048;
leny=200*0.3048;
lenz=50*0.3048;%m
dx1=10*0.3048;
dy1=10*0.3048;
dz=lenz;%m
time_total=2*365*24*60*60;%sec
dt=1*24*60*60;%sec
Nt=round(time_total/dt);

Nx=lenx/dx1;
Ny=leny/dy1;
dx=zeros(Ny,Nx);
dy=zeros(Ny,Nx);
dx(:,:)=dx1;
dy(:,:)=dx1;

%  Nx=20;
%  Ny=20;

% dx1=randi([10,20],[Ny,Nx]);
% dy1=randi([10,20],[Ny,Nx]);
%  dx=dx1.*0.3048;
%  dy=dy1.*0.3048;

rw=0.3*0.3048;%meter
re=sqrt(dx.*dy/pi);%m
p1 =-4.167;
p2=12.5;
p3=-15.46;
p4=7.125;

pc_f=@(s) (p1*s.^3 + p2*s.^2 + p3*s+p4)*6894.7573;%pa
dpc_dsw=@(s) (3.*p1*s.^2 + 2.*p2*s + p3)*6894.7573;

perm=randi([10,350],[Ny,Nx]);
% perm=xlsread('perm.xlsx');
perm=perm(1:Ny,1:Nx)*0.001*9.869233e-13;%m2
% perm=ones(Ny,Nx)*0.5*9.869233e-13;
bo0=1;
bw0=1;

swi=0.2;
swc=0.2;
sor=0.2;

Cw=1/6894.7573*1e-5;
Co=1/6894.7573*1e-5;
Cr=1/6894.7573*1e-8;

fi=xlsread('fi.xlsx');
fi=fi(1:Ny,1:Nx);

Bo=@(p) 1./(1+Co*(p-3200*6894.7573));
Bw=@(p) 1./(1+Cw*(p-3200*6894.7573));

muob=3.5;
pbo=2635;
mu_o=@(p) (muob*(p/pbo)^(2.6*(p^1.187)*exp(-11.513+(-8.98e-5*p))))/10^3;

muwb=3.5;
pbw=2635;
mu_w=@(p) (muob*(p/pbw)^(2.6*(p^1.187)*exp(-11.513+(-8.98e-5*p))))/10^3;

init_p=5000*6894.7573;%pa
pressure=zeros(Ny,Nx,Nt);
sat=zeros(Ny,Nx,Nt);
pressure(:,:,1)=init_p;
sat(:,:,1)=swi;


zarayeb=zeros(Nx*Ny,Nx*Ny);
sabet=zeros(Nx*Ny,1);


trans=@(kr1,kr2,mu1,mu2,B1,B2,pre1,pre2,k1,k2,d1,d2) 2*(kr1/mu1/B1*(pre1>=pre2)+kr2/mu2/B2*(pre2>pre1))/(d1)/(d1/k1+d2/k2);

q_inj=15*5.615*(0.3048)^3/86400;
pbh=3600*6894.7573;



for n=2:Nt
    n %#ok

    for j=1:Ny
    
        for i=1:Nx
            index=i+(j-1)*Nx;
            landa_o=rel_perm(sat(j,i,n-1),2)/mu_o(pressure(j,i,n-1))/Bo(pressure(j,i,n-1));
            landa_w=rel_perm(sat(j,i,n-1),1)/mu_w(pressure(j,i,n-1))/Bw(pressure(j,i,n-1));
            cpooi=fi(j,i)*(1-sat(j,i,n-1))/dt*(Cr/Bo(pressure(j,i,n-1))+Co);
            cswoi=-fi(j,i)/dt/Bo(pressure(j,i,n-1));
            cpowi=fi(j,i)*(sat(j,i,n-1))/dt*(Cr/Bw(pressure(j,i,n-1))+Cw);
            cswwi=fi(j,i)/dt/Bw(pressure(j,i,n-1))-dpc_dsw(sat(j,i,n-1))*cpowi;
            alpha=-cswwi/cswoi;
            
            
            if i==1&&j==1
            
                trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
%                 trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_r+trans_y_o_d+cpooi)-alpha*(trans_x_w_r+trans_y_w_d+cpowi);
                B=trans_x_o_r+alpha*trans_x_w_r;
%                 C=trans_x_o_l+alpha*trans_x_w_l;
%                 D=trans_y_o_u+alpha*trans_y_w_u;
                E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_r*(pc_f(sat(j,i+1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_d*(pc_f(sat(j+1,i,n-1))-pc_f(sat(j,i,n-1))));
                F=F-alpha*abs(q_inj)/dx(j,i)/dy(j,i)/dz;
                zarayeb(index,index+1)=B;
%                 zarayeb(index,index-1)=C;
%                 zarayeb(index,index-Nx)=D;
                zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;
 
                
                
             elseif i==1 && j==Ny
                
                trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
%                 trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
%                 trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_r+trans_y_o_u+cpooi)-alpha*(trans_x_w_r+trans_y_w_u+cpowi);
                B=trans_x_o_r+alpha*trans_x_w_r;
%                 C=trans_x_o_l+alpha*trans_x_w_l;
                D=trans_y_o_u+alpha*trans_y_w_u;
%                 E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_r*(pc_f(sat(j,i+1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_u*(pc_f(sat(j-1,i,n-1))-pc_f(sat(j,i,n-1))));
                zarayeb(index,index+1)=B;
%                 zarayeb(index,index-1)=C;
                zarayeb(index,index-Nx)=D;
%                 zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;


            elseif i==Nx && j==1
                
%                 trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
%                 trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_l+trans_y_o_d+cpooi)-alpha*(trans_x_w_l+trans_y_w_d+cpowi);
%                 B=trans_x_o_r+alpha*trans_x_w_r;
                C=trans_x_o_l+alpha*trans_x_w_l;
%                 D=trans_y_o_u+alpha*trans_y_w_u;
                E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_l*(pc_f(sat(j,i-1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_d*(pc_f(sat(j+1,i,n-1))-pc_f(sat(j,i,n-1))));
%                 zarayeb(index,index+1)=B;
                zarayeb(index,index-1)=C;
%                 zarayeb(index,index-Nx)=D;
                zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;

            elseif i==Nx && j==Ny
                
%                 trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
%                 trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
%                 trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_l+trans_y_o_u+cpooi)-alpha*(trans_x_w_l+trans_y_w_u+cpowi);
%                 B=trans_x_o_r+alpha*trans_x_w_r;
                C=trans_x_o_l+alpha*trans_x_w_l;
                D=trans_y_o_u+alpha*trans_y_w_u;
%                 E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_l*(pc_f(sat(j,i-1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_u*(pc_f(sat(j-1,i,n-1))-pc_f(sat(j,i,n-1))));
                WC=2*pi*perm(j,i)*dz/log(re(j,i)/rw);            
                A=A-WC/dx(j,i)/dy(j,i)/dz*landa_o-alpha*WC/dx(j,i)/dy(j,i)/dz*landa_w;
                F=F-WC/dx(j,i)/dy(j,i)/dz*landa_o*pbh-alpha*WC/dx(j,i)/dy(j,i)/dz*landa_w*pbh;
%                 zarayeb(index,index+1)=B;
                zarayeb(index,index-1)=C;
                zarayeb(index,index-Nx)=D;
%                 zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;

                


            elseif i==1
                
                trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_r+trans_y_o_u+trans_y_o_d+cpooi)-alpha*(trans_x_w_r+trans_y_w_u+trans_y_w_d+cpowi);
                B=trans_x_o_r+alpha*trans_x_w_r;
%                 C=trans_x_o_l+alpha*trans_x_w_l;
                D=trans_y_o_u+alpha*trans_y_w_u;
                E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_r*(pc_f(sat(j,i+1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_u*(pc_f(sat(j-1,i,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_d*(pc_f(sat(j+1,i,n-1))-pc_f(sat(j,i,n-1))));
                zarayeb(index,index+1)=B;
%                 zarayeb(index,index-1)=C;
                zarayeb(index,index-Nx)=D;
                zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;

            
            elseif i==Nx
     
%                 trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_l+trans_y_o_u+trans_y_o_d+cpooi)-alpha*(trans_x_w_l+trans_y_w_u+trans_y_w_d+cpowi);
%                 B=trans_x_o_r+alpha*trans_x_w_r;
                C=trans_x_o_l+alpha*trans_x_w_l;
                D=trans_y_o_u+alpha*trans_y_w_u;
                E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_l*(pc_f(sat(j,i-1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_u*(pc_f(sat(j-1,i,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_d*(pc_f(sat(j+1,i,n-1))-pc_f(sat(j,i,n-1))));
%                 zarayeb(index,index+1)=B;
                zarayeb(index,index-1)=C;
                zarayeb(index,index-Nx)=D;
                zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;

            elseif j==1
                
                trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
%                 trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_r+trans_x_o_l+trans_y_o_d+cpooi)-alpha*(trans_x_w_r+trans_x_w_l+trans_y_w_d+cpowi);
                B=trans_x_o_r+alpha*trans_x_w_r;
                C=trans_x_o_l+alpha*trans_x_w_l;
%                 D=trans_y_o_u+alpha*trans_y_w_u;
                E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_r*(pc_f(sat(j,i+1,n-1))-pc_f(sat(j,i,n-1)))+trans_x_w_l*(pc_f(sat(j,i-1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_d*(pc_f(sat(j+1,i,n-1))-pc_f(sat(j,i,n-1))));
                zarayeb(index,index+1)=B;
                zarayeb(index,index-1)=C;
%                 zarayeb(index,index-Nx)=D;
                zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;

            elseif j==Ny
                
                trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
%                 trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
%                 trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_r+trans_x_o_l+trans_y_o_u+cpooi)-alpha*(trans_x_w_r+trans_x_w_l+trans_y_w_u+cpowi);
                B=trans_x_o_r+alpha*trans_x_w_r;
                C=trans_x_o_l+alpha*trans_x_w_l;
                D=trans_y_o_u+alpha*trans_y_w_u;
%                 E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_r*(pc_f(sat(j,i+1,n-1))-pc_f(sat(j,i,n-1)))+trans_x_w_l*(pc_f(sat(j,i-1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_u*(pc_f(sat(j-1,i,n-1))-pc_f(sat(j,i,n-1))));
                zarayeb(index,index+1)=B;
                zarayeb(index,index-1)=C;
                zarayeb(index,index-Nx)=D;
%                 zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;

            
            else
                trans_x_o_r=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i+1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i+1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i+1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i+1,n-1)),pressure(j,i,n-1),pressure(j,i+1,n-1),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_o_l=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j,i-1,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j,i-1,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j,i-1,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j,i-1,n-1)),pressure(j,i,n-1),pressure(j,i-1,n-1),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_o_u=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j-1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j-1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j-1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j-1,i,n-1)),pressure(j,i,n-1),pressure(j-1,i,n-1),perm(j,i),perm(j-1,i),dx(j,i),dy(j-1,i));
                trans_y_o_d=trans(rel_perm(sat(j,i,n-1),2),rel_perm(sat(j+1,i,n-1),2),mu_o(pressure(j,i,n-1)),mu_o(pressure(j+1,i,n-1)),Bo(pressure(j,i,n-1)),Bo(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n-1)),mu_w(pressure(j+1,i,n-1)),Bw(pressure(j,i,n-1)),Bw(pressure(j+1,i,n-1)),pressure(j,i,n-1),pressure(j+1,i,n-1),perm(j,i),perm(j+1,i),dx(j,i),dy(j+1,i));
                A=-(trans_x_o_r+trans_x_o_l+trans_y_o_u+trans_y_o_d+cpooi)-alpha*(trans_x_w_r+trans_x_w_l+trans_y_w_u+trans_y_w_d+cpowi);
                B=trans_x_o_r+alpha*trans_x_w_r;
                C=trans_x_o_l+alpha*trans_x_w_l;
                D=trans_y_o_u+alpha*trans_y_w_u;
                E=trans_y_o_d+alpha*trans_y_w_d;
                F=-(cpooi+alpha*cpowi)*pressure(j,i,n-1)+alpha*(trans_x_w_r*(pc_f(sat(j,i+1,n-1))-pc_f(sat(j,i,n-1)))+trans_x_w_l*(pc_f(sat(j,i-1,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_u*(pc_f(sat(j-1,i,n-1))-pc_f(sat(j,i,n-1)))+trans_y_w_d*(pc_f(sat(j+1,i,n-1))-pc_f(sat(j,i,n-1))));
                zarayeb(index,index+1)=B;
                zarayeb(index,index-1)=C;
                zarayeb(index,index-Nx)=D;
                zarayeb(index,index+Nx)=E;
                zarayeb(index,index)=A;
                sabet(index,1)=F;
            end
        end
    end
    
    javab=linsolve(zarayeb,sabet);
    
    for j=1:Ny
        for i=1:Nx
            index=i+(j-1)*Nx;
            pressure(j,i,n)=javab(index);
        end
    end
    
    
    for j=1:Ny
        for i=1:Nx
            if i==1 && j==1
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
                sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_r*(pressure(j,i+1,n)-pc_f(sat(j,i+1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_d*(pressure(j+1,i,n)-pc_f(sat(j+1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));

                sat(j,i,n)=sat(j,i,n)+1/cswwi/dx(j,i)/dy(j,i)/dz*q_inj; 
            elseif  i==1 && j==Ny
                
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                 trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
%                 trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
              sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_r*(pressure(j,i+1,n)-pc_f(sat(j,i+1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_u*(pressure(j-1,i,n)-pc_f(sat(j-1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
            elseif i==Nx && j==1
%                 trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
            sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_l*(pressure(j,i-1,n)-pc_f(sat(j,i-1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_d*(pressure(j+1,i,n)-pc_f(sat(j+1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
            elseif i==Nx && j==Ny
              
%                 trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
%                 trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
              sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_l*(pressure(j,i-1,n)-pc_f(sat(j,i-1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_u*(pressure(j-1,i,n)-pc_f(sat(j-1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
               WC=2*pi*perm(j,i)*dz/log(re(j,i)/rw);            
               sat(j,i,n)=sat(j,i,n)-1/cswwi*WC/dx(j,i)/dy(j,i)/dz*landa_w*(pressure(j,i,n)-pbh);
            
            elseif i==1
    
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
%                 trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
               sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_r*(pressure(j,i+1,n)-pc_f(sat(j,i+1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_u*(pressure(j-1,i,n)-pc_f(sat(j-1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_d*(pressure(j+1,i,n)-pc_f(sat(j+1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
            elseif i==Nx
%                 trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
             sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_l*(pressure(j,i-1,n)-pc_f(sat(j,i-1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_u*(pressure(j-1,i,n)-pc_f(sat(j-1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_d*(pressure(j+1,i,n)-pc_f(sat(j+1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
            elseif j==1
               
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
%                 trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
               sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_r*(pressure(j,i+1,n)-pc_f(sat(j,i+1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_x_w_l*(pressure(j,i-1,n)-pc_f(sat(j,i-1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_d*(pressure(j+1,i,n)-pc_f(sat(j+1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
            elseif j==Ny
             
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
%                 trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
             sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_r*(pressure(j,i+1,n)-pc_f(sat(j,i+1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_x_w_l*(pressure(j,i-1,n)-pc_f(sat(j,i-1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_u*(pressure(j-1,i,n)-pc_f(sat(j-1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
            else
                trans_x_w_r=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i+1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i+1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i+1,n)),pressure(j,i,n),pressure(j,i+1,n),perm(j,i),perm(j,i+1),dx(j,i),dx(j,i+1));
                trans_x_w_l=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j,i-1,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j,i-1,n)),Bw(pressure(j,i,n)),Bw(pressure(j,i-1,n)),pressure(j,i,n),pressure(j,i-1,n),perm(j,i),perm(j,i-1),dx(j,i),dx(j,i-1));
                trans_y_w_u=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j-1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j-1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j-1,i,n)),pressure(j,i,n),pressure(j-1,i,n),perm(j,i),perm(j-1,i),dx(j,i),dx(j-1,i));
                trans_y_w_d=trans(rel_perm(sat(j,i,n-1),1),rel_perm(sat(j+1,i,n-1),1),mu_w(pressure(j,i,n)),mu_w(pressure(j+1,i,n)),Bw(pressure(j,i,n)),Bw(pressure(j+1,i,n)),pressure(j,i,n),pressure(j+1,i,n),perm(j,i),perm(j+1,i),dx(j,i),dx(j+1,i));
                sat(j,i,n)=sat(j,i,n-1)+1/cswwi*(trans_x_w_r*(pressure(j,i+1,n)-pc_f(sat(j,i+1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_x_w_l*(pressure(j,i-1,n)-pc_f(sat(j,i-1,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_u*(pressure(j-1,i,n)-pc_f(sat(j-1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))+trans_y_w_d*(pressure(j+1,i,n)-pc_f(sat(j+1,i,n-1))-(pressure(j,i,n)-pc_f(sat(j,i,n-1))))-cpowi*(pressure(j,i,n)-pressure(j,i,n-1)));
            end
        end
    end
end



for n=1:Nt
    pcolor(real(sat(:,:,n)));
    drawnow
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1