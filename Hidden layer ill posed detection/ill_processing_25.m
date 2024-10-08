clc
clear

load('PLOT4.mat')
IW_30 = IWW(1:25,:);
B_30 = BB(1:25,:);
BiasMatrix_30 = repmat(B_30,1,Q);
tempH_30 = IW_30 * Pn_train + BiasMatrix_30;
H_30 = 1 ./ (1 + exp(-tempH_30));
A_30=cond(H_30);
LW_30 =pinv(H_30') * Tn_train';
t_30 = elmpredict(Pn_test,IW_30,B_30,LW_30,'sig',0);
T_30=mapminmax('reverse',t_30,ou);
L= length(P_test);
r_30=(L*sum(T_30.*T_test)-sum(T_30)*sum(T_test))/sqrt(((L*sum((T_30).^2)-(sum(T_30))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
rmse_30=sqrt(sum((T_test-T_30).^2)/L);
%% Tikhonov
[U_L2,D_L2,V_L2]=csvd(H_30');
reg_param=l_curve(U_L2,D_L2,Tn_train');
RR=eye(25);
A_L2=cond(H_30*H_30'+(reg_param)*RR);
LW_L2=inv(H_30*H_30'+(reg_param)*RR)*H_30*Tn_train';
t_L2 = elmpredict(Pn_test,IW_30,B_30,LW_L2,'sig',0);
T_L2=mapminmax('reverse',t_L2,ou);
L= length(P_test);
r_L2=(L*sum(T_L2.*T_test)-sum(T_L2)*sum(T_test))/sqrt(((L*sum((T_L2).^2)-(sum(T_L2))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
rmse_L2=sqrt(sum((T_test-T_L2).^2)/L);
%% TSVD
[U_SVD,S_SVD,V_SVD]=svd(H_30');
s_SVD=diag(pinv(S_SVD));
LW_SVD=zeros(size(LW_30,1),1);
rmse_SVD=[];
rmse_TSVD=[];
      for i=1:10
          LWi_SVD=V_SVD(:,i)*s_SVD(i)*U_SVD(:,i)'*Tn_train';
          LW_SVD=LWi_SVD+LW_SVD;
          %% training set
          t_SVD = elmpredict(Pn_train,IW_30,B_30,LW_SVD,'sig',0);
          T_SVD = mapminmax('reverse',t_SVD,ou);
          rmse_S=sqrt(sum((T_train-T_SVD).^2)/60);
            rmse_SVD=[rmse_SVD;rmse_S];
            %% test set
            t_TSVD = elmpredict(Pn_test,IW_30,B_30,LW_SVD,'sig',0);
            T_TSVD = mapminmax('reverse',t_TSVD,ou);
             rmse_TS=sqrt(sum((T_test-T_TSVD).^2)/L);      
            rmse_TSVD=[rmse_TSVD;rmse_TS];
      end    
      [C,D]=min(rmse_TSVD);
      D=1;
      rmse_TSVD(D);
      A_TSVD=S_SVD(1,1)/S_SVD(D,D);      
%% whitening transformation
c=cov(H_30'); %求矩阵x的协方差矩阵
[F,V]=eigs(c); %F的列向量为对应特征向量,V为最大特征值对角阵
v=V^(-0.5);
W=(F*v*(F)');
H_WT=W*H_30;
A_WT=cond(H_WT);
L_WT =pinv(H_WT') * Tn_train'; 
BiasMatrix_WT = repmat(B_30,1,22);
tempH_WT = IW_30 * Pn_test + BiasMatrix_WT;
H_WT2 =W*( 1 ./ (1 + exp(-tempH_WT)));
t_WT = (H_WT2' * L_WT)';
T_WT = mapminmax('reverse',t_WT,ou);
r_30=(L*sum(T_30.*T_test)-sum(T_30)*sum(T_test))/sqrt(((L*sum((T_30).^2)-(sum(T_30))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
rmse_WT=sqrt(sum((T_test-T_WT).^2)/L);
% Optionally, if the number of nodes is high, it can be combined with TSVD
% [U,S,V]=svd(H_WT');
% s=diag(pinv(S));
% LW_WT=zeros(size(LW_30,1),1);
% rmse_WT=[];
%       for i=1:6
%           LWi=V(:,i)*s(i)*U(:,i)'*Tn_train';
%           LW_WT=LWi+LW_WT;
%           BiasMatrix_WT = repmat(B_30,1,22);
%           tempH_WT = IW_30 * Pn_test + BiasMatrix_WT;
%           H_WT2 =W*( 1 ./ (1 + exp(-tempH_WT)));
%           t_WT = (H_WT2' * LW_WT)';
%           T_WT = mapminmax('reverse',t_WT,ou);
%           rmse_W=sqrt(sum((T_test-T_WT).^2)/L);
%             rmse_WT=[rmse_WT;rmse_W];
%       end
% min(rmse_WT)
% S(1,1)/S(2,2)
