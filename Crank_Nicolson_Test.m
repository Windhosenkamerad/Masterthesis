clc;
clear all;
close all;


R = 100;

q = 0.005;
k = 100;
dt = 10;
dr = 10;

N = R/dr;

a = k*dt/dr;

nOnes = ones(N, 1); %Hilfe für Diagonalelemente    

T1 = 170;

%initiales Temperaturprofil
for r = dr:dr:R
    T(r/dr,1) = 150;
end

for t=dt:dt:10000           %Timeloop
    for r = dr:dr:R - dr    %spatial loop
    
a = k*dt/(2*r*dr);

T(1,1) = T1;

alpha = a - (a*r/dr);   %Faktoren der Tridiagonalmatrix
beta = 1 + 2*a*r/dr;
gamma = a + (a*r/dr);

%Erstellen der Tridiagonalmatrix
A = diag(beta * nOnes, 0) + diag(alpha*nOnes(1:N-1), -1) - diag(gamma*nOnes(1:N-1), 1);

d(1,1) = T(1,t/dt);                             %Dirichlet Quellvektor
d(N,1) = T(N-2,t/dt) - T(N,t/dt) + 4*q*dr/k;    %Neumann Quellvektor

%Erstellen des Quellvektors
if r/dr > 1
d(r/dr,1) = -alpha*T(r/dr - 1,t/dt) + (2-beta)*T(r/dr,t/dt) + gamma*T(r/dr + 1,t/dt);

A(1,1) = 1;         %Dirichlet
A(1,2) = 0;         %Dirichlet
A(N,N-2) = 1;
A(N,N-1) = 0;       %Neumann
A(N,N) = -1;        %Neumann

end
                        T(:,t/dt+1) = inv(A) * d;

    end
                    
end

disp('calculations finished');

plot(1:dr:R,T(:,1));
hold on
plot(1:dr:R,T(:,end))

%contourf(T)
