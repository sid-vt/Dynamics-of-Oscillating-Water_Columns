%% Initializing code
clear;clc;close all;

%% Defining System Parameters
% mass, tension, length, alpha and beta
m = 5;
T = 15;
L = 0.1;
k = 10;

% For proportional damping
alpha = 0.010;
beta = 0.050;

L1 = L;
L2 = L;
L3 = L;
L4 = L;
L5 = L;
L6 = L;
L7 = L;

m1 = m;
m2 = m;
m3 = m;
m4 = m;
m5 = m;
m6 = m;

k1 = k;
k2 = k;
k3 = k;
k4 = k;
k5 = k;
k6 = k;

%% Defining Mass, Stiffness, Damping Matrix
M = [m1, 0, 0, 0, 0, 0; 0, m2, 0, 0, 0, 0; 0, 0, m3, 0, 0, 0; 0, 0,...
   0, m4, 0, 0; 0, 0, 0, 0, m5, 0; 0, 0, 0, 0, 0, m6];

K = [k1 + T/L1 + T/L2, -(T/L2), 0, 0, 0, 0; -(T/L2),...
  k2 + T/L2 + T/L3, -(T/L3), 0, 0, 0; 0, -(T/L3),...
  k3 + T/L3 + T/L4, -(T/L4), 0, 0; 0, 0, -(T/L4),...
  k4 + T/L4 + T/L5, -(T/L5), 0; 0, 0, 0, -(T/L5),...
  k5 + T/L5 + T/L6, -(T/L6); 0, 0, 0, 0, -(T/L6), k6 + T/L6 + T/L7];

C = alpha*M + beta*K;

%% Eigen Analysis

[Evec, Eval] = eig(K,M);
%Sorting required
[~,index] = sort(diag(Eval));

%% Results of Eigen Analysis:
% lambda, U(mass normalized eigen vector)
lambda = diag(Eval(index,index));
U = Evec(:,index);

%% Mode Shapes
dof = 6;
masscnt = linspace(1,dof,dof);

figure('Name','Mode Shape','NumberTitle','off')
for i=1:dof
    subplot(dof,1,i)
    plot(masscnt,U(:,i)','-or','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
    yline(0,'-.k','LineWidth',1.5);
    yticks(-2:0.5:2)
    xlabel('x');
    ylabel('H');
    grid on;
end

%% Data required
% Natural frequencies from sqrt(lambda):
wn = sqrt(lambda);

% Zeta by considering proportional Damping:
zeta = diag(U'*C*U)./wn/2;

% Damped frequency:
wd = sqrt(1-zeta.^2).*wn;

%% Conditions and Forces
Q0 = [0;0;0];   % No forces
w = 0;
x0=[10,10,10,10,10,10];
v0=[0,0,0,0,0,0];
t = 0:0.01:5;

%% Free Vibration

for i=1:dof
    xt(i,:) = exp(-zeta(i)*wn(i)*t).*(x0(i)*cos(wd(i)*t)+...
        (v0(i)+zeta(i)*wn(i)*x0(i))/wd(i)*sin(wd(i)*t));
    
    figure(2)
    plot(t,xt(i,:));hold on;

end
figure (2)
title('Free Frequency Response for individual mass')
legend('Mass 1','Mass 2','Mass 3','Mass 4','Mass 5','Mass 6')
grid on;

xt = U(:,1)*xt(1,:) + U(:,2)*xt(2,:) + U(:,3)*xt(3,:) +...
    U(:,4)*xt(4,:) + U(:,5)*xt(5,:) + U(:,6)*xt(6,:);
    
figure('Name','Free Frequency combined Response','NumberTitle','off')
plot(t,xt)
title('Free Frequency combined Response')
legend('Mass 1','Mass 2','Mass 3','Mass 4','Mass 5','Mass 6')
grid on;

