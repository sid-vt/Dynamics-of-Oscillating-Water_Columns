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
% masscnt = linspace(1,dof,dof);
% figure('Name','Mode Shape','NumberTitle','off')
% for i=1:dof
%     subplot(dof,1,i)
%     plot(masscnt,U(:,i)','-or','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
%     yline(0,'-.k','LineWidth',1.5);
%     yticks(-2:0.5:2)
%     xlabel('x');
%     ylabel('H');
%     grid on;
% end

%% Data required
% Natural frequencies from sqrt(lambda):
wn = sqrt(lambda);

% Zeta by considering proportional Damping:
zeta = diag(U'*C*U)./wn/2;

% Damped frequency:
wd = sqrt(1-zeta.^2).*wn;

%%
% Excitation Frequency
w = 0.5;
t = 0:0.01:100;
Q0 = [1;2;3;4;5;6];

%% From Mathematica computation 
etacos1 = (exp(-t*wn(1)*zeta(1)).*(-wd(1)*(-w^2 + wd(1)^2 + wn(1)^2*zeta(1)^2)*cos(t*wd(1)) +... 
   exp(t*wn(1)*zeta(1))*wd(1).*((-w^2 + wd(1)^2 + wn(1)^2*zeta(1)^2)*cos(t*w) +... 
      2*w*wn(1)*zeta(1)*sin(t*w)) - wn(1)*zeta(1)*(w^2 + wd(1)^2 + wn(1)^2*zeta(1)^2)*sin(t*wd(1))))/(((w - wd(1))^2 +... 
   wn(1)^2*zeta(1)^2)*((w + wd(1))^2 + wn(1)^2*zeta(1)^2));

etasin2 = (exp(-t*wn(2)*zeta(2)).*(2*w*wd(2)*wn(2)*zeta(2)*cos(t*wd(2)) +...
    exp(t*wn(2)*zeta(2))*wd(2).*(-2*w*wn(2)*zeta(2)*cos(t*w) +...
    (-w^2 + wd(2)^2 + wn(2)^2*zeta(2)^2)*sin(t*w)) + w*(w^2 - wd(2)^2 +...
    wn(2)^2*zeta(2)^2)*sin(t*wd(2))))/(((w - wd(2))^2 + wn(2)^2*zeta(2)^2)*((w + wd(2))^2 + wn(2)^2*zeta(2)^2));

etacos3 = (exp(-t*wn(3)*zeta(3)).*(-wd(3)*(-w^2 + wd(3)^2 + wn(3)^2*zeta(3)^2)*cos(t*wd(3)) +... 
   exp(t*wn(3)*zeta(3))*wd(3).*((-w^2 + wd(3)^2 + wn(3)^2*zeta(3)^2)*cos(t*w) +... 
      2*w*wn(3)*zeta(3)*sin(t*w)) - wn(3)*zeta(3)*(w^2 + wd(3)^2 + wn(3)^2*zeta(3)^2)*sin(t*wd(3))))/(((w - wd(3))^2 +... 
   wn(3)^2*zeta(3)^2)*((w + wd(3))^2 + wn(3)^2*zeta(3)^2));

etasin4 = (exp(-t*wn(4)*zeta(4)).*(2*w*wd(4)*wn(4)*zeta(4)*cos(t*wd(4)) +...
    exp(t*wn(4)*zeta(4))*wd(4).*(-2*w*wn(4)*zeta(4)*cos(t*w) +...
    (-w^2 + wd(4)^2 + wn(4)^2*zeta(4)^2)*sin(t*w)) + w*(w^2 - wd(4)^2 +...
    wn(4)^2*zeta(4)^2)*sin(t*wd(4))))/(((w - wd(4))^2 + wn(4)^2*zeta(4)^2)*((w + wd(4))^2 + wn(4)^2*zeta(4)^2));

etacos5 = (exp(-t*wn(5)*zeta(5)).*(-wd(5)*(-w^2 + wd(5)^2 + wn(5)^2*zeta(5)^2)*cos(t*wd(5)) +... 
   exp(t*wn(5)*zeta(5))*wd(5).*((-w^2 + wd(5)^2 + wn(5)^2*zeta(5)^2)*cos(t*w) +... 
      2*w*wn(5)*zeta(5)*sin(t*w)) - wn(5)*zeta(5)*(w^2 + wd(5)^2 + wn(5)^2*zeta(5)^2)*sin(t*wd(5))))/(((w - wd(5))^2 +... 
   wn(5)^2*zeta(5)^2)*((w + wd(5))^2 + wn(5)^2*zeta(5)^2));

etasin6 = (exp(-t*wn(6)*zeta(6)).*(2*w*wd(6)*wn(6)*zeta(6)*cos(t*wd(6)) +...
    exp(t*wn(6)*zeta(6))*wd(6).*(-2*w*wn(6)*zeta(6)*cos(t*w) +...
    (-w^2 + wd(6)^2 + wn(6)^2*zeta(6)^2)*sin(t*w)) + w*(w^2 - wd(6)^2 +...
    wn(6)^2*zeta(6)^2)*sin(t*wd(6))))/(((w - wd(6))^2 + wn(6)^2*zeta(6)^2)*((w + wd(6))^2 + wn(6)^2*zeta(6)^2));


eta1 = ((U(:,1)'*Q0)/wd(1))*etacos1;
eta2 = ((U(:,2)'*Q0)/wd(2))*etasin2;
eta3 = ((U(:,3)'*Q0)/wd(3))*etacos3;
eta4 = ((U(:,1)'*Q0)/wd(1))*etasin4;
eta5 = ((U(:,2)'*Q0)/wd(2))*etacos5;
eta6 = ((U(:,3)'*Q0)/wd(3))*etasin6;
q = eta1.*U(:,1)+eta2.*U(:,2)+eta3.*U(:,3)+eta4.*U(:,4)+eta5.*U(:,5)+eta6.*U(:,6);

figure(2)
plot(t,q)
legend('Mass 1','Mass 2','Mass 3','Mass 4','Mass 5','Mass 6')
grid on;