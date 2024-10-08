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
N=10;
IW = rand(N,R) * 2 - 1;
B = rand(N,1);
BiasMatrix = repmat(B,1,Q);
tempH = IW * Pn_train + BiasMatrix;
H = 1 ./ (1 + exp(-tempH));
A=cond(H);
LW0 =pinv(H') * Tn_train';
t = elmpredict(Pn_test,IW,B,LW0,'sig',0);
t_xunlian = elmpredict(Pn_train,IW,B,LW0,'sig',0);
T=mapminmax('reverse',t,ou);
L= length(P_test);
r=(L*sum(T.*T_test)-sum(T)*sum(T_test))/sqrt(((L*sum((T).^2)-(sum(T))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
rmse=sqrt(sum((T_test-T).^2)/L);

%% Setting the gain threshold
threshold=0;
%% Calculate the original network loss
t_xunlian = elmpredict(Pn_train,IW,B,LW0,'sig',0);
T_xunlian=mapminmax('reverse',t_xunlian,ou);
L_xunlian= length(T_xunlian);
original_loss=sum((T_train-T_xunlian).^2);
pruned_indices = [];
num_nodes = N;
    for node_index = 1:num_nodes
        % Remove a node
        pruned_IW = IW;
        pruned_B = B;
        pruned_IW(node_index,:) = [];
        pruned_B(node_index,:) = [];
        pruned_BiasMatrix = repmat(pruned_B,1,Q);
        pruned_tempH = pruned_IW * Pn_train + pruned_BiasMatrix;
        pruned_H = 1 ./ (1 + exp(-pruned_tempH));
        pruned_LW =pinv(pruned_H') * Tn_train';
        % Calculating network loss after pruning
        pruned_t = elmpredict(Pn_train,pruned_IW,pruned_B,pruned_LW,'sig',0);
        pruned_T=mapminmax('reverse',pruned_t,ou);
        pruned_loss=sum((T_train-pruned_T).^2);
        % Calculate the information gain, if the original loss is higher, this node will be removed
        information_gain =  pruned_loss-original_loss;
        % If the information gain is below the threshold, pruning
        if information_gain < threshold
            pruned_indices = [pruned_indices, node_index];
        end
    end
       pruned_IW = IW;
        pruned_B = B;
        pruned_IW(pruned_indices,:) = [];
        pruned_B(pruned_indices,:) = [];
        pruned_BiasMatrix = repmat(pruned_B,1,Q);
         pruned_tempH = pruned_IW * Pn_train + pruned_BiasMatrix;
        pruned_H = 1 ./ (1 + exp(-pruned_tempH));
        pruned_LW =pinv(pruned_H') * Tn_train';
        pruned_t = elmpredict(Pn_test,pruned_IW,pruned_B,pruned_LW,'sig',0);
        pruned_T=mapminmax('reverse',pruned_t,ou);
        L= length(P_test);
        pruned_r=(L*sum(pruned_T.*T_test)-sum(pruned_T)*sum(T_test))/sqrt(((L*sum((pruned_T).^2)-(sum(pruned_T))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)))
        pruned_rmse=sqrt(sum((T_test-pruned_T).^2)/L)
        MAPE1=mean(abs((T_test - pruned_T) ./ T_test)) * 100;
        MAE1=mean(abs(T_test - pruned_T))
        
save('GP-ELM.mat')