%create a low rank signal
l=100; %spectral bands
n=1000; % pixels
shots=20; %shots
rank=8; %rank of the signal

A=randn(l,rank);
B=rand(n,rank);
F=A*B'; %input signal

%calculate and substract mean
f=mean(F,2);
F1=F-(kron(f,ones(1,n)));

%%Create sensing matrices

% random sensing matrix
Q=rand(shots,100);
% designed matrix
transmittance = 1/1; %maximun possible transmittance
[Q_tilde ,Sigmarest]= designbinary_pca(shots,F1*F1'./n,transmittance);

% sense and reconstruct using Pseudo-Inverse

Yr=Q*F; %random measurements
Yd=Q_tilde*F; %designed measurments

Fr=pinv(Q)*Yr;
Fd=pinv(Q_tilde)*Yd;

%% comparison
pix=round(n*rand()); %pixel to be plotted
plot(1:l,F(:,pix),'red',1:l,Fr(:,pix),'blue',1:l,Fd(:,pix),'black'),
legend('Ground-truth','Random','Designed')

fprintf('||F-Fd||_2=%f, ||F-Fr||_2=%f, SNR_Designed=%f, SNR_Random=%f\n',norm(F-Fd),norm(F-Fr),snr(F,F-Fd),snr(F,F-Fr));

%% comparison of covariance

showSigmas=[F1*F1'./n,zeros(l,2)-1,Sigmarest,zeros(l,2)-1,F1*F1'./n-Sigmarest];

figure,imagesc(showSigmas),xlabel('Original Covariance, Information Deflated, Residual')

