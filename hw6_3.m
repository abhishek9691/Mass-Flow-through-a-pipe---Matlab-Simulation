function OneD_transient_mass_transfer
 
close all
 
D = 2E-9; %m^2/s, OG, Pe = 0.1
%D = 2E-10; %m^2/s, OG, Pe = 1
%D = 2E-11; %m^2/s, OG, Pe = 10
%D = 2E-12; %m^2/s, OG, Pe = 100

V = 2E-8;%m/s
Cin = 0; %mM
nw = 0; %mol/s*m^2
%L = 0.01; %m
%L = 0.1; %Pe = 1
L = 1; %Pe = 10
%L = 10; %Pe = 100
diam = 0.001; %m
k = 1;
P = 2*3.14*diam*L; 
A = 2*3.14*(diam/2)*L + 2*3.14*(diam/2)^2;
Pe = (V*L)/D 

Nx=201;    %number of points in x direction
Nt=101;    %number of points in time
m=0;       %cartesian
%gamma=(nw*P*(L^2))/(Cin*D*A);
gamma = 0; 
psi_w=1;
psi_0=1;
 
Xmax=1;
Tmax=3;
dX=Xmax/(Nx-1);
dT=Tmax/(Nt-1);
 
options=[];
T=linspace(0,Tmax,Nt);
X=linspace(0,Xmax,Nx);
pde=@(x,t,C,DCDx)pdepb(x,t,C,DCDx,gamma);
IC=@(x)pbIC(x,psi_0);
BC=@(yl,Cl,yr,Cr,t)pbBC(yl,Cl,yr,Cr,t,psi_w);
sol=pdepe(m,pde,IC,BC,X,T,options);  %solution      
C=sol(:,:,1);                                         %extract solution 
 
x_value1=2;
x_value2=81;   
x_value3=121; 
x_value4=181;
plot(T,C(:,x_value1),T,C(:,x_value2),T,C(:,x_value3),T,C(:,x_value4))   %Plot vs t at specified x locations
xlabel('T')
ylabel('Concentration')
xlim([0 0.75]);
legend(['X = ' num2str(x_value1*dX-dX)],[num2str(x_value2*dX-dX)],[num2str(x_value3*dX-dX)],[num2str(x_value4*dX-dX)])
 
figure(2)   %Plot vs x at specified times
t1=1;
t2=11;
t3=41;
t4=101;
plot(X,C(t1,:),X,C(t2,:),X,C(t3,:),X,C(t4,:))
xlabel('X')
ylabel('Concentration')
%ylim([0 1.5])
%xlim([0 Xmax])
legend(['T = ' num2str(t1*dT-dT)],[num2str(t2*dT-dT)],[num2str(t3*dT-dT)], [num2str(t4*dT-dT)])

%dimensional form 
x = X.*L;
C_dimension = C.*Cin; 
t = (T.*(L^2))/D; 

figure(3)
plot(t,C_dimension(:,x_value1),t,C_dimension(:,x_value2),t,C_dimension(:,x_value3),t,C_dimension(:,x_value4))   %Plot vs t at specified x locations
xlabel('t')
ylabel('Concentration')
legend(['X = ' num2str(x_value1*dX-dX)],[num2str(x_value2*dX-dX)],[num2str(x_value3*dX-dX)],[num2str(x_value4*dX-dX)])

figure(4)
plot(x,C_dimension(t1,:),x,C_dimension(t2,:),x,C_dimension(t3,:),x,C_dimension(t4,:))
xlabel('x')
ylabel('Concentration')
%ylim([0 1.5])
%xlim([0 L])
legend(['T = ' num2str(t1*dT-dT)],[num2str(t2*dT-dT)],[num2str(t3*dT-dT)], [num2str(t4*dT-dT)])

function [c,f,s]=pdepb(x,t,C,DCDx,gamma)  %define pde system – see Matlab help for pdepe
c=1;           %coefficient on time derivative
f=DCDx;  
s=gamma;     
 
function C0=pbIC(x,psi_0)   %initial condition
C0=0.1*x;   %mol/cm^3
 
function [pl,ql,pr,qr]=pbBC(yl,Cl,yr,Cr,t,psi_w)  %boundary conditions
pl=0;   %left side bc. This one sets the flux to zero
ql=1;       
pr=Cr-psi_w;      %right side bc. This one sets the concentration to psi_w
qr=0;
