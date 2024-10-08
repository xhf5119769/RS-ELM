clc
clear
load('data.mat')
load('y.mat')
data=data/100;
X=data';
X_mean=mean(X);
QQ = randperm(size(X,1));
P_train = X(QQ(1:60),:);
T_train = y(QQ(1:60),:);
P_test = X(QQ(61:end),:);
T_test = (y(QQ(61:end),:));
[Z,W,P,T] = dosc(P_train,T_train,3,1E-4);
[Ztest,Ttest] = dosc_pred(P_test,W,P);
X1=Z';
X2=Ztest';

X=[X1 X2];
=[T_train;T_test];