%create a low rank signal
l=100;
n=1000;
shots=12;

A=randn(l,shots);
B=rand(n,shots);
F=A*B';
%calculate and substract mean
f=mean(F,2);
F1=F-(kron(f,ones(1,n)));

%%Create sensing matrices

% random sensing matrix
Q=rand(shots,100);
% designed matrix
Q_tilde = designbinary_pca(shots,l,F1*F1'./n);

% sense and reconstruct using Pseudo-Inverse

Yr=Q*F;
Yd=Q_tilde*F;

Fr=pinv(Q)*Yr;
Fd=pinv(Q_tilde)*Yd;

%% comparison
pix=560; %pixel to be plotted
plot(1:l,F(:,pix),'red',1:l,Fr(:,pix),'blue',1:l,Fd(:,pix),'black'),
legend('Ground-truth','Random','Designed')

