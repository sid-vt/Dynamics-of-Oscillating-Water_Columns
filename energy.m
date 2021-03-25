%% Initializing code
clear;clc;close all;

%% Defining System Parameters
% mass, tension, length, alpha and beta
m = 200;
T = 50;
L = 2;
k = 04;

% For proportional damping
alpha = 0.10;
beta = 0.50;

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

% Excitation Frequency
w = 0.5;
t = 0:0.01:100;
Q0 = [1;2;3;4;5;6];

%%
G1 = 1./sqrt((1-(w/wn(1))^2)^2 + (2*zeta(1)*(w/wn(1)))^2);
G2 = 1./sqrt((1-(w/wn(2))^2)^2 + (2*zeta(2)*(w/wn(2)))^2);
G3 = 1./sqrt((1-(w/wn(3))^2)^2 + (2*zeta(3)*(w/wn(3)))^2);
G4 = 1./sqrt((1-(w/wn(4))^2)^2 + (2*zeta(4)*(w/wn(4)))^2);
G5 = 1./sqrt((1-(w/wn(5))^2)^2 + (2*zeta(5)*(w/wn(5)))^2);
G6 = 1./sqrt((1-(w/wn(6))^2)^2 + (2*zeta(6)*(w/wn(6)))^2);

phi1 = atan((2*zeta(1)*(w/wn(1)))/(1-(w/wn(1))));
phi2 = atan((2*zeta(2)*(w/wn(2)))/(1-(w/wn(2))));
phi3 = atan((2*zeta(3)*(w/wn(3)))/(1-(w/wn(3))));
phi4 = atan((2*zeta(4)*(w/wn(4)))/(1-(w/wn(4))));
phi5 = atan((2*zeta(5)*(w/wn(5)))/(1-(w/wn(5))));
phi6 = atan((2*zeta(6)*(w/wn(6)))/(1-(w/wn(6))));

%%
q1 = (U(:,1)'*Q0/(wn(1))^2)*G1*U(:,1)*exp(1i*(w*t-phi1));
q2 = (U(:,2)'*Q0/(wn(2))^2)*G2*U(:,2)*exp(1i*(w*t-phi2));
q3 = (U(:,3)'*Q0/(wn(3))^2)*G3*U(:,3)*exp(1i*(w*t-phi3));
q4 = (U(:,4)'*Q0/(wn(4))^2)*G4*U(:,4)*exp(1i*(w*t-phi4));
q5 = (U(:,5)'*Q0/(wn(5))^2)*G5*U(:,5)*exp(1i*(w*t-phi5));
q6 = (U(:,6)'*Q0/(wn(6))^2)*G6*U(:,6)*exp(1i*(w*t-phi6));

q = real(q1)+imag(q2)+real(q3)+imag(q4)+real(q5)+imag(q6);

qdot1=((G1*w*U(:,1)'*Q0/wn(1)^2)/(sqrt((1 - w^2/wn(1)^2)^2 + (4*w^2*zeta(1)^2)/wn(1)^2)))*U(:,1)*1i*exp(1i*(-phi1 + t*w));
qdot2=((G2*w*U(:,2)'*Q0/wn(2)^2)/(sqrt((1 - w^2/wn(2)^2)^2 + (4*w^2*zeta(2)^2)/wn(2)^2)))*U(:,1)*1i*exp(2i*(-phi2 + t*w));
qdot3=((G3*w*U(:,3)'*Q0/wn(3)^2)/(sqrt((1 - w^2/wn(3)^2)^2 + (4*w^2*zeta(3)^2)/wn(3)^2)))*U(:,1)*1i*exp(1i*(-phi3 + t*w));
qdot4=((G4*w*U(:,4)'*Q0/wn(4)^2)/(sqrt((1 - w^2/wn(4)^2)^2 + (4*w^2*zeta(4)^2)/wn(4)^2)))*U(:,1)*1i*exp(1i*(-phi4 + t*w));
qdot5=((G5*w*U(:,5)'*Q0/wn(5)^2)/(sqrt((1 - w^2/wn(5)^2)^2 + (4*w^2*zeta(5)^2)/wn(5)^2)))*U(:,1)*1i*exp(1i*(-phi5 + t*w));
qdot6=((G6*w*U(:,6)'*Q0/wn(6)^2)/(sqrt((1 - w^2/wn(6)^2)^2 + (4*w^2*zeta(6)^2)/wn(6)^2)))*U(:,1)*1i*exp(1i*(-phi6 + t*w));

qdot = real(qdot1)+imag(qdot2)+real(qdot3)+imag(qdot4)+real(qdot5)+imag(qdot6);

for i=1:length(t)
    energy1(:,i) = 0.5*K*q1(:,i).^2+0.5*M*qdot1(:,i).^2;
    energy2(:,i) = 0.5*K*q2(:,i).^2+0.5*M*qdot2(:,i).^2;
    energy3(:,i) = 0.5*K*q3(:,i).^2+0.5*M*qdot3(:,i).^2;
    energy4(:,i) = 0.5*K*q4(:,i).^2+0.5*M*qdot4(:,i).^2;
    energy5(:,i) = 0.5*K*q5(:,i).^2+0.5*M*qdot5(:,i).^2;
    energy6(:,i) = 0.5*K*q6(:,i).^2+0.5*M*qdot6(:,i).^2;
end

for i=1:length(t)
    energycombine(:,i) = 0.5*K*q(:,i).^2+0.5*M*qdot(:,i).^2;
end

figure(2)
plot(t,energy6)

figure(3)
plot(t,energycombine(6,:))
