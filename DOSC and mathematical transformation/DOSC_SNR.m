clc
clear
load('data.mat')
load('y.mat')
data=data/100;
X=data';
X_mean=mean(X);
SNR=[];
DATA_X={};
for i=1:11
% X_Z=X-ones(82,1)*X_mean;
[Z,W,P,T] = dosc(X,y,i,1E-4);
X_DOSC=Z';
DATA_X{i}=X_DOSC;
noise=data-X_DOSC;
DATA=reshape(X_DOSC,984,1);
NO=reshape(noise,984,1);
snr_value = snr(DATA, NO);
SNR(i)=snr_value;
end
plot(SNR)
ylabel('SNR (dB)');
xlabel('Number of orthogonal components');

