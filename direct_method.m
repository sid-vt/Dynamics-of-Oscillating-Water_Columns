%% Initializing code
clear;clc;close all;

%% Defining System Parameters
% mass, tension, length, alpha and beta
m1 = 5;m2 = 5;m3 = 5;m4 = 5;m5 = 5;
T = 15;
L = 10;

%% Defining Mass, Stiffness, Damping Matrix

M = diag([m1,m2,m3]);
K=[10 -6 0; -6 14 -8;0 -8 8];

% proportional damping
alpha = 0.010;
beta = 00.050;
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
dof = 3;
masscnt = linspace(1,dof,dof);

figure('Name','Mode Shape','NumberTitle','off')
for i=1:dof
    subplot(dof,1,i)
    plot(masscnt,U(:,i)','-or','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
    yline(0,'-.k','LineWidth',1.5);
    ylim([-2,2])
    yticks(-2:0.5:2)
    xlabel('x');
    ylabel('H');
    grid on;
end

%Draw

%% Data required
% Natural frequencies from sqrt(lambda):
wn = sqrt(lambda);

% Zeta by considering proportional Damping:
zeta = diag(U'*C*U)./wn/2;

% Damped frequency:
wd = sqrt(1-zeta.^2).*wn;

w = 1;

Q0 = [1;0;0];
t = 0:0.1:100;
%%

G1 = 1./sqrt((1-(w/wd(1))^2)^2 + (2*zeta(1)*(w/wd(1)))^2);
G2 = 1./sqrt((1-(w/wd(2))^2)^2 + (2*zeta(2)*(w/wd(2)))^2);
G3 = 1./sqrt((1-(w/wd(3))^2)^2 + (2*zeta(3)*(w/wd(3)))^2);

phi1 = atan((2*zeta(1)*(w/wd(1)))/(1-(w/wd(1))));
phi2 = atan((2*zeta(2)*(w/wd(2)))/(1-(w/wd(2))));
phi3 = atan((2*zeta(3)*(w/wd(3)))/(1-(w/wd(3))));
% 
% q1 = (U(:,1)'*Q0/(wd(1))^2)*G1*U(:,1)*exp(1i*(w*t-phi1));
% q2 = (U(:,2)'*Q0/(wd(2))^2)*G2*U(:,2)*exp(1i*(w*t-phi2));
% q3 = (U(:,3)'*Q0/(wd(3))^2)*G3*U(:,3)*exp(1i*(w*t-phi3));
% 
% qdot1 = G1*U(:,1)'*Q0*w/wn(1)^2*U(:,1)*(1i*exp(1i*(t*w - phi1)));
% qdot2 = G2*U(:,2)'*Q0*w/wn(2)^2*U(:,2)*(1i*exp(1i*(t*w - phi2)));
% qdot3 = G3*U(:,3)'*Q0*w/wn(3)^2*U(:,3)*(1i*exp(1i*(t*w - phi3)));
% 
% q = real(q1)+real(q2)+real(q3);
% qdot = real(qdot1)+real(qdot2)+real(qdot3);

q1 = (U(:,1)'*Q0/(wd(1))^2)*G1*U(:,1)*cos(w*t-phi1);
q2 = (U(:,2)'*Q0/(wd(2))^2)*G2*U(:,2)*cos(w*t-phi2);
q3 = (U(:,3)'*Q0/(wd(3))^2)*G3*U(:,3)*cos(w*t-phi3);

qdot1 = G1*U(:,1)'*Q0*w/wn(1)^2*U(:,1)*(w*(sin(t*w - phi1)));
qdot2 = G2*U(:,2)'*Q0*w/wn(2)^2*U(:,2)*(w*(sin(t*w - phi2)));
qdot3 = G3*U(:,3)'*Q0*w/wn(3)^2*U(:,3)*(w*(sin(t*w - phi3)));

q = q1+q2+q3;
qdot = qdot1+qdot2+qdot3;

figure(2)
plot(t,q,'LineWidth',1.5)
legend('1','2','3')

w = linspace(0,5,1000);
for i=1:length(w)
    G11(i) = 1./sqrt((1-(w(i)/wd(1))^2)^2 + (2*zeta(1)*(w(i)/wd(1)))^2);
    G22(i) = 1./sqrt((1-(w(i)/wd(2))^2)^2 + (2*zeta(2)*(w(i)/wd(2)))^2);
    G33(i) = 1./sqrt((1-(w(i)/wd(3))^2)^2 + (2*zeta(3)*(w(i)/wd(3)))^2);
end
figure (3)
plot(w,G11);hold on;
plot(w,G22);hold on;
plot(w,G33);

for i=1:length(t)
    energy(:,i) = 0.5*(K*q(:,i).^2)+0.5*(M*qdot(:,i).^2);
end

figure (4)
plot(t,energy(1,:));hold on;
plot(t,energy(2,:));hold on;
plot(t,energy(3,:));

for i=1:length(t)
    energy(:,i) = 0.5*(q(:,i)'*K*q(:,i).^2)+0.5*(qdot(:,i)'*M*qdot(:,i).^2);
end

figure (5)
plot(t,energy(1,:));hold on;
plot(t,energy(2,:));hold on;
plot(t,energy(3,:));
