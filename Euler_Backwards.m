clc;
clear all;
close all;

yr = 60*60*24*365.5;
R = 10000;
time = 10000*yr;

k = 50;
rho = 5000;
cp = 800;

q1 =  -0.00005;
q2 = 0.5;
kappa = k/(rho*cp);
dt = 100*yr;
dr = 100;

%allocate space for matrices
T = zeros(R/dr, time/dt);
d = zeros(R/dr,1);

N = R/dr;

nOnes = ones(N, 1); %Hilfe für Diagonalelemente    

T1 = 170;           % Temperatur an unterer Schichtgrenze -> Dirichlet RB

%initiales Temperaturprofil
for r = dr:dr:R
    %T(r/dr,1) = T1 - r/dr; %mit küsntlichem Gradient
    T(r/dr,1) = T1 + 0.1*r;
end
%T(1,1) = T1;

for t = dt:dt:time           %Timeloop
    for r = dr:dr:R        %spatial loop
    
a = kappa*dt/(r*dr);

alpha = a*r/dr - a;   %Faktoren der Tridiagonalmatrix
beta = 1 + 2*a*r/dr;
gamma = a*r/dr + a;

%% Erstellen der Tridiagonalmatrix
A = diag(beta * nOnes, 0) - diag(alpha*nOnes(1:N-1), -1) - diag(gamma*nOnes(1:N-1), 1);

%% Quellvektor mit Neumann & Dirichlet RB
d(r/dr,1) = T(r/dr,t/dt);         %Erstellen des Quellvektors
%d(1,1) = T(1,t/dt) + alpha*T1;   %Dirichlet
d(1,1) = T1;                      %Dirichlet
d(N,1) = T(N,t/dt) - 2*dr*q1/k;    %Neumann

A(1,1) = 1;                       %Dirichlet
A(1,2) = 0;                       %Dirichlet
A(N,N-1) = - alpha - gamma;       %Neumann
A(N,N) = beta;                    %Neumann

%% Einstellungen für Neumann-Bedingungen an beiden Enden der Schicht
% A(1,1) = beta;                    
% A(1,2) = -alpha - gamma;                  
% A(N,N-1) = -alpha - gamma;
% A(N,N) = beta;
% 
% d(r/dr,1) = T(r/dr,t/dt);
% d(1,1) = T(1,t/dt) + alpha*2*dr*q1/k;
% d(N,1) = T(N,t/dt) + gamma*2*dr*q2/k;

                        T(:,t/dt+1) = inv(A) * d;

    end
                    
end

disp('calculations finished');

%% Testplots
figure(1)
contourf(T)
colorbar
ylabel('Radius in dr')
xlabel('Zeit in dt')

figure(2)
plot(T(:,1));
figure(3)
plot((T(:,end)));