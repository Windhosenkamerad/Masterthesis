clc;
clear all;
close all;

N = 20;
n = 50;
dx = 1;

r = 1:1:N;
q = 5;
k = 100;
dt = 10;
dr = 10;

a = k*dt/dr;
% alpha = a/dr;
% beta = 1 - 2*a/dr - 2*a/r;
% gamma = 2*a/r + a/dr;

nOnes = ones(N, 1) ;    %Erschaffen einer Tridiagonalmatrix und Korrektur des letzten Eintrags
%A = diag(beta * nOnes, 0) + diag(alpha*nOnes(1:n-1), -1) + diag(gamma*nOnes(1:n-1), 1);
%B = (1/dx^2)*diag(4 * nOnes, 0) + diag(nOnes(1:n-1), -1) + diag(nOnes(1:n-1), 1);


T1 = 4200;

%initiales Temperaturprofil
for r = 1:1:N
    T(r,1) = 4000;
    %T(r,2) = 4000 - r*70;
end

for t=dt:dt:10000    %Timeloop
    for r = dr:dr:1000   %spatial loop
    
a = k*dt/(2*r*dr);    
alpha = a - (a*r/dr);             %Faktoren der Tridiagonalmatrix
beta = 2*a*r/dr +1;
gamma = a+(a*r/dr);
j = a*r/dr - a;
k = 1 - 2*a*r/dr;
l = a + a*r/dr;

A = diag(beta * nOnes, 0) + diag(alpha*nOnes(1:N-1), -1) - diag(gamma*nOnes(1:N-1), 1);
B = diag(k * nOnes, 0) + diag(j*nOnes(1:N-1), -1) + diag(l*nOnes(1:N-1), 1);

%Einbauen der Dirichlet/Neumann-RB in der Koeffizientenmatrix
A(N,N) = 1;
A(N, N-1) = -1;
A(1,1) = 1;
A(1,2) = 0;
T(1,1) = T1;       % T(r1) = T1 !
T(N,1) = -q*dr/k;  % T(N,t) - T(N-1,t) = q*dt --> letzte Gleichung der Matrix-Vektormultiplikation
B(1,2) = 0;
B(N,N-1) = 0;

           
            BT = B * T(:,t/dt);

            T(:,t/dt +1 ) = inv(A) * BT;
    end
end

disp('calculations finished');

% figure(1)
% plot(1:size(T,1),T)
% title('old')
% 
% figure(2)
% title('new')
% plot(1:size(Tnew,1),Tnew)


contourf(T)
