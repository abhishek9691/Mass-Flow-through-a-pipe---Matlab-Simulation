numx = 101;   %number of grid points in x
numt = 2000;  %number of time steps to be iterated over 
dx = 1/(numx - 1);
dt = 0.00005;

x = 0:dx:1;   %vector of x values, to be used for plotting

C = zeros(numx,numt);   %initialize everything to zero

%specify initial conditions
t(1) = 0;      %t=0
C(1,1) = 0;    %C=0 at x=0
C(1,numx) = 0; %C=0 at x=1
mu = 0.5;
sigma = 0.05;
for i=2:numx-1
   C(i,1) = exp(-(x(i)-mu)^2/(2*sigma^2)) / sqrt(2*pi*sigma^2);
end

%iterate difference equation - note that C(1,j) and C(numx,j) always remain 0
for j=1:numt
   t(j+1) = t(j) + dt;
   for i=2:numx-1
      C(i,j+1) = C(i,j) + (dt/dx^2)*(C(i+1,j) - 2*C(i,j) + C(i-1,j)); 
   end
end
close all
figure(1);
hold on;
plot(x,C(:,1));
%plot(x,C(:,11));
plot(x,C(:,101));
%plot(x,C(:,1001));
plot(x,C(:,2001));
xlabel('x');
xlim([0.5 1])
ylabel('c(x,t)');

