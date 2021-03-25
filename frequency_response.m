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
% 
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

w = linspace(0,15,1000);
for i=1:length(w)
    G1(i) = 1./sqrt((1-(w(i)/wd(1))^2)^2 + (2*zeta(1)*(w(i)/wd(1)))^2);
    G2(i) = 1./sqrt((1-(w(i)/wd(2))^2)^2 + (2*zeta(2)*(w(i)/wd(2)))^2);
    G3(i) = 1./sqrt((1-(w(i)/wd(3))^2)^2 + (2*zeta(3)*(w(i)/wd(3)))^2);
    G4(i) = 1./sqrt((1-(w(i)/wd(4))^2)^2 + (2*zeta(4)*(w(i)/wd(4)))^2);
    G5(i) = 1./sqrt((1-(w(i)/wd(5))^2)^2 + (2*zeta(5)*(w(i)/wd(5)))^2);
    G6(i) = 1./sqrt((1-(w(i)/wd(6))^2)^2 + (2*zeta(6)*(w(i)/wd(6)))^2);
    
    phi1(i) = atan((2*zeta(1)*(w(i)/wd(1)))/(1-(w(i)/wd(1))));
    phi2(i) = atan((2*zeta(2)*(w(i)/wd(2)))/(1-(w(i)/wd(2))));
    phi3(i) = atan((2*zeta(3)*(w(i)/wd(3)))/(1-(w(i)/wd(3))));
    phi4(i) = atan((2*zeta(4)*(w(i)/wd(4)))/(1-(w(i)/wd(4))));
    phi5(i) = atan((2*zeta(5)*(w(i)/wd(5)))/(1-(w(i)/wd(5))));
    phi6(i) = atan((2*zeta(6)*(w(i)/wd(6)))/(1-(w(i)/wd(6))));
end

figure (3)
plot(w,G1);hold on;
plot(w,G2);hold on;
plot(w,G3);hold on;
plot(w,G4);hold on;
plot(w,G5);hold on;
plot(w,G6);
grid on;
legend('Mass 1','Mass 2','Mass 3','Mass 4','Mass 5','Mass 6')

figure (4)
plot(w,phi1);hold on;
plot(w,phi2);hold on;
plot(w,phi3);hold on;
plot(w,phi4);hold on;
plot(w,phi5);hold on;
plot(w,phi6);
grid on;
legend('Mass 1','Mass 2','Mass 3','Mass 4','Mass 5','Mass 6')
