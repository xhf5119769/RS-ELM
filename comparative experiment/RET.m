clc
clear
load('data.mat')
load('y.mat')
data=data/100;
y=y';
QQ = randperm(size(data,2));
 load('QQ.mat')
P_train = data(:,QQ(1:60));
T_train = y(:,QQ(1:60));
P_test = data(:,QQ(61:end));
T_test = y(QQ(:,61:end));
%% DOSC 
[Z,W,P,TT] = dosc(P_train',T_train',4,1E-4);
[Ztest,Ttest] = dosc_pred(P_test',W,P);
Z=Z';
Ztest=Ztest';
%% First-order derivative spectra
window_size = 3;  % ���ڴ�С��������
order = 1;  % ����
sg_coeff = sgolay(order, window_size);
for i=1:60
RFD = conv(Z(:,i), sg_coeff(:, order+1), 'same');
P_train_FOD(:,i)=RFD;
end

for i=1:22
RFD = conv(Ztest(:,i), sg_coeff(:, order+1), 'same');
P_test_FOD(:,i)=RFD;
end
%% min-max normalization
[Pn_train,in] = mapminmax(P_train_FOD);
[Tn_train,ou] = mapminmax(T_train);
Pn_test = mapminmax('apply',P_test_FOD,in);

M    =  50; % Number of Extra_Trees in the ensemble
k    =  3;  % Number of attributes selected to perform the random splits 
            % 1 <k <= total number of attributes 
nmin =  10;  % Minimum number of points for each leaf
[ensemble,TT_train] = buildAnEnsemble(M,k,nmin,[Pn_train',Tn_train'],0);
% Run the ensemble on a validation dataset
TT = predictWithAnEnsemble(ensemble,[Pn_test',T_test'],0);
T = mapminmax('reverse',TT',ou);
N=length(T);
RMSE_RBF=sqrt(sum((T_test-T).^2)/N)
R2=sqrt((N*sum(T.*T_test)-sum(T)*sum(T_test))/sqrt(((N*sum((T).^2)-(sum(T))^2)*(N*sum((T_test).^2)-(sum(T_test))^2))))
MAE1=mean(abs(T_test - T));
MAPE1=mean(abs((T_test - T) ./ T_test)) * 100;


