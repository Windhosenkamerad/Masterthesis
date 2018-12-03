function [ Tnew ] = WLG(T, r_cmb, dt, dr, q, T1 ,r1 , rhog, cp, k)

%Erstellen des Quellvektors
d = T;

N = size(T,1); %N=size(T,1)+1
%dr = (r_cmb - r1)/N;      %Einteilen der Schicht in N Gridpunkte
%dr = size(T,1);

kappa = k/(rhog*cp);

nOnes = ones(N, 1) ;      %Erschaffen einer Tridiagonalmatrix und Korrektur des letzten Eintrags

    for r = r1:dr:r_cmb   %spatial loop
    
a = kappa*dt/(r*dr);

alpha = a*r/dr - a;   %Faktoren der Tridiagonalmatrix
beta = 1 + 2*a*r/dr;
gamma = a*r/dr + a;

%Erstellen der Tridiagonalmatrix
A = diag(beta * nOnes, 0) - diag(alpha*nOnes(1:N-1), -1) - diag(gamma*nOnes(1:N-1), 1);


%d(1,1) = T(1,t/dt) + alpha*T1;
d(1,1) = T1(end);
%d(N,1) = T(N-1,1) - gamma*2*dr*q/k; %old
d(N,1) = T(N-1,1) - 2*dr*q/k; %new

A(1,1) = 1;
A(1,2) = 0;
A(N,N-1) = - alpha - gamma;        %Neumann
A(N,N) = beta;                     %Neumann

                        Tnew = A\d;
                        
                        if Tnew(Tnew<0)
                            disp('Warning: negative temperatures!')
                        end

    end
end


