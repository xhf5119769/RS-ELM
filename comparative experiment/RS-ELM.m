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
N=20;
IW = rand(N,R) * 2 - 1;
B = rand(N,1);
BiasMatrix = repmat(B,1,Q);
tempH = IW * Pn_train + BiasMatrix;
H = 1 ./ (1 + exp(-tempH));
A=cond(H);
LW0 =pinv(H') * Tn_train';
t = elmpredict(Pn_test,IW,B,LW0,'sig',0);
T=mapminmax('reverse',t,ou);
L= length(P_test);
r1=(L*sum(T.*T_test)-sum(T)*sum(T_test))/sqrt(((L*sum((T).^2)-(sum(T))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)))
rmse=sqrt(sum((T_test-T).^2)/L)
%%  H is discretized
[m,n]=size(H);
 H_max=max(max(H));
H_min=min(min(H));

C=(H_max-H_min)/9;
U=[1:1:10]';
Z=[];
for i=1:m
    for j=1:n
        Z(i,j)=floor((H(i,j)-H_min)/C)+1;
    end
end

%% Turning discretized H into an information system
BB=zeros(m,n);
for ii=1:m
    for jj=1:n
        BB(ii,Z(ii,jj))=1;
    end
end
BB=BB(:,1:10)';

r=[];
r_mean=[];
R_mean=[];
%% Calculate the dependency degree
i=1;j=1;
for i=1:N
    for j=1:N
r(i,j)=DependencyDegree(U,BB(:,i),BB(:,j));
    end
end

i =1; j=1;

%% Dependency degree averaging
for i=1:N
    for j=1:N
        r_mean(i,j)=(r(i,j)+r(j,i))./2;
    end
end
R_mean=tril(r_mean,-1);
%% 
IW2=IW;
B2=B;
%% Delete nodes according to rules
 % Setting Thresholds
 threshold=0.35;
% Iteration matrix
i=1;
j=1;
for j = 1:size(R_mean, 1)
    for i = i+1:size(R_mean, 2)
        % First check if (i, j) is greater than the threshold value
        if R_mean(i, j) > threshold
            % Rule 1
           if nnz(R_mean(i, :) == 0) == numel(R_mean(i, :)) - 1 & nnz(R_mean(:, j) == 0) == numel(R_mean(:, j)) - 1
            % Delete a node at random
            IW2(i,:)=0;
            B2(i,:)=0;
           end 
            %Rule 2
            count1 = sum(R_mean(i, :) > threshold);
            count2 = sum(R_mean(:, j) > threshold);
            if count1>1 & count2 >1
              %Both nodes are deleted
             IW2(i,:)=0;
            B2(i,:)=0;
             IW2(j,:)=0;
            B2(j,:)=0;            
            end 
               %Rule 3
            if count1>1 & count2 < 2 
                %Delete the first node
            IW2(i,:)=0;
            B2(i,:)=0;     
            end
            if count1<2 & count2 >1 
                %Delete the second node
            IW2(j,:)=0;
            B2(j,:)=0;   
            end
            %Rule 4
            if count1<2 & count2 <2
                if sum(R_mean(i, :)) > sum(R_mean(:, j));
            IW2(i,:)=0;
            B2(i,:)=0;   
                else
             IW2(j,:)=0;
            B2(j,:)=0;     
                end
            end
        end
    end
end

% Deleting rows with all zeros using logical indexes
nonzero_rows = any(IW2 ~= 0, 2);
IW2 = IW2(nonzero_rows, :);
nonzero_rows = any(B2 ~= 0, 2);
B2 = B2(nonzero_rows, :);

BiasMatrix2 = repmat(B2,1,Q);

tempH2 = IW2 * Pn_train + BiasMatrix2;
H2 = 1 ./ (1 + exp(-tempH2));
TF='sig';
TYPE=0;
LW2 =pinv(H2') * Tn_train';
t = elmpredict(Pn_test,IW2,B2,LW2,TF,TYPE);
T=mapminmax('reverse',t,ou);
r2=(L*sum(T.*T_test)-sum(T)*sum(T_test))/sqrt(((L*sum((T).^2)-(sum(T))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
rmse2=sqrt(sum((T_test-T).^2)/L);
r2
r1

NN=size(B2,1);
%% WT and TSVD
c=cov(H2'); 
[F,V]=eigs(c); 
v=V^(-0.5);
W=(F*v*(F)');
H_WT=W*H2;
L_WT =pinv(H_WT') * Tn_train'; 
BiasMatrix_WT = repmat(B2,1,22);
tempH_WT = IW2 * Pn_test + BiasMatrix_WT;
H_WT2 =W*( 1 ./ (1 + exp(-tempH_WT)));
t_WT = (H_WT2' * L_WT)';
T_WT = mapminmax('reverse',t_WT,ou);
r=(L*sum(T_WT.*T_test)-sum(T_WT)*sum(T_test))/sqrt(((L*sum((T_WT).^2)-(sum(T_WT))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)))
rmse_WT=sqrt(sum((T_test-T_WT).^2)/L)

% i=1;
% [U,S,V]=svd(H_WT');
% s=diag(pinv(S));
% LW_WT=zeros(size(LW2,1),1);
% rmse_WT=[];
%       for i=1:NN
%           LWi=V(:,i)*s(i)*U(:,i)'*Tn_train';
%           LW_WT=LWi+LW_WT;
%           BiasMatrix_WT = repmat(B2,1,22);
%           tempH_WT = IW2 * Pn_test + BiasMatrix_WT;
%           H_WT2 =W*( 1 ./ (1 + exp(-tempH_WT)));
%           t_WT = (H_WT2' * LW_WT)';
%           T_WT = mapminmax('reverse',t_WT,ou);
%           rmse_W=sqrt(sum((T_test-T_WT).^2)/L);
%           rmse_WT=[rmse_WT;rmse_W];
%           MAE1(i)=mean(abs(T_test - T_WT));
%          MAPE1(i)=mean(abs((T_test - T_WT) ./ T_test)) * 100;
%           r_WT(i)=sqrt((L*sum(T_WT.*T_test)-sum(T_WT)*sum(T_test))/sqrt(((L*sum((T_WT).^2)-(sum(T_WT))^2)*(L*sum((T_test).^2)-(sum(T_test))^2))));
%       end
% r3=max(r_WT)
% rmse_f=min(rmse_WT);
% MAE_F=min(MAE1);
% MAPE_F=min(MAPE1)

save('RS-ELM.mat')