clc
clear
load('data.mat')
load('y.mat')
data=data/100;
X=data';
X_mean=mean(X);
% X_Z=X-ones(82,1)*X_mean;
[Z,W,P,T] = dosc(X,y,4,1E-4);
X_DOSC=Z';
save('X_DOSC',"X_DOSC")