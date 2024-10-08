clc
clear
load('data.mat')
load('y.mat')
data=data/100;
y=y';
QQ = randperm(size(data,2));
 % load('QQ.mat')
P_train = data(:,QQ(1:60));
T_train = y(:,QQ(1:60));
P_test = data(:,QQ(61:end));
T_test = y(QQ(:,61:end));
%% DOSC 
[Z,W,P,T] = dosc(P_train',T_train',4,1E-4);
[Ztest,Ttest] = dosc_pred(P_test',W,P);
Z=Z';
Ztest=Ztest';
%% First-order derivative spectra
window_size = 3;  % 窗口大小（奇数）
order = 1;  % 阶数
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
%% Modeling of the ELM
[R,Q]=size(Pn_train);
N=5;
IW = rand(N,R) * 2 - 1;
B = rand(N,1);
BiasMatrix = repmat(B,1,Q);
tempH = IW * Pn_train + BiasMatrix;
H = 1 ./ (1 + exp(-tempH));
A=cond(H);
[U_L2,D_L2,V_L2]=csvd(H');
reg_param=l_curve(U_L2,D_L2,Tn_train');
RR=eye(N);
LW0=inv(H*H'+reg_param*RR)*H*Tn_train';
t = elmpredict(Pn_test,IW,B,LW0,'sig',0);
t_xunlian = elmpredict(Pn_train,IW,B,LW0,'sig',0);
T=mapminmax('reverse',t,ou);
L= length(P_test);
r=(L*sum(T.*T_test)-sum(T)*sum(T_test))/sqrt(((L*sum((T).^2)-(sum(T))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
rmse=sqrt(sum((T_test-T).^2)/L);
MAE=mean(abs(T_test - T));
MAPE=mean(abs((T_test - T) ./ T_test)) * 100;
%% 
IWW=IW;
BB=B;
HH=H;
LWW=LW0;
RMSE=[];
R2=[];
CD=[];
R2=[R2;r];
RMSE=[RMSE;rmse];
CD=[CD;A];
%% loop, adding one node at a time
for i=1:40
    IW_new=rand(1,R) * 2 - 1;
    B_new = rand(1,1);
    IWW=[IWW;IW_new];
    BB=[BB;B_new];
    BiasMatrix_new = repmat(B_new,1,Q);
    tempH_new = IW_new * Pn_train + BiasMatrix_new;
    H_new = 1 ./ (1 + exp(-tempH_new));
    HH=[HH;H_new];
    A_new=cond(HH);
    LW_1 =pinv(H_new') * (Tn_train-t_xunlian)';
    LWW=[LWW;LW_1];
    t_new = elmpredict(Pn_test,IWW,BB,LWW,'sig',0);
    t_xunlian = elmpredict(Pn_train,IWW,BB,LWW,'sig',0);
    T_new=mapminmax('reverse',t_new,ou);
    r_new=sqrt((L*sum(T_new.*T_test)-sum(T_new)*sum(T_test))/sqrt(((L*sum((T_new).^2)-(sum(T_new))^2)*(L*sum((T_test).^2)-(sum(T_test))^2))));
    rmse_new=sqrt(sum((T_test-T_new).^2)/L);
    MAE(i+1)=mean(abs(T_test - T_new));
    MAPE(i+1)=mean(abs((T_test - T_new) ./ T_test)) * 100;
    R2=[R2;r_new];
    RMSE=[RMSE;rmse_new];
    CD=[CD;A_new];
end
[rmse_f,D]=min(RMSE)
MAE_F=MAE(D)
MAPE_F=MAPE(D)
R2_F=R2(D)

save('IR-ELM.mat')