%calculator

%N=number of grid=indexb
%ju=j of top
%jd=j of bottom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
clc;
N=10000;
Nx=100;
Ny=100;
jd=floor((N-1)/Nx)+1;
ju=Ny-jd+1;
i=N-(jd-1)*Nx;

disp(['ju=',num2str(ju)]);
disp(['i=',num2str(i)]);