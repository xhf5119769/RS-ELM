clc
clear
load('data.mat')
load('y.mat')
data=data/100;
X=data';
X_mean=mean(X);
E_finl=[];
for A=1:5
QQ = randperm(size(X,1));
P_train = X(QQ(1:60),:);
T_train = y(QQ(1:60),:);
P_test = X(QQ(61:end),:);
T_test = (y(QQ(61:end),:))';
EE=[];
for i=0:11
[Z,W,P,T] = dosc(P_train,T_train,i,1E-4);
[Ztest,Ttest] = dosc_pred(P_test,W,P);
Z=Z';
Ztest=Ztest';
T_t=T_train';
[Pn_train,in] = mapminmax(Z);
[Tn_train,ou] = mapminmax(T_t);
Pn_test = mapminmax('apply',Ztest,in);
E=[];
for j=1:5
[IW,B,LW,TF,TYPE] = elmtrain(Pn_train,Tn_train,20);
Tn_test = elmpredict(Pn_test,IW,B,LW,TF,TYPE);
%∑¥πÈ“ªªØ
T = mapminmax('reverse',Tn_test,ou);
N=length(T);
e=sqrt(sum((T_test-T).^2)/N);
E=[E,e];
R2=sqrt((N*sum(T.*T_test)-sum(T)*sum(T_test))^2/((N*sum((T).^2)-(sum(T))^2)*(N*sum((T_test).^2)-(sum(T_test))^2)))
end
E_mean=mean(E);
EE=[EE,E_mean];
end
E_finl=[E_finl;EE];
end
E_F=mean(E_finl);

plot(E_F(1:9))
ylabel('RMSE of the test set');
xlabel('Number of orthogonal components');