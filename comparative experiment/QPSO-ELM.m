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


%% parameter setting
m=40; 
q=50; 
xi=20;
sigma=1;
beta=0.1;
[R,Q]=size(Pn_train);

for i=1:m
 pop(i,:)=round(40*rand(1));
 IW{i,:} = rand(pop(i),R) * 2 - 1;
B{i,:}= rand(pop(i),1);
BiasMatrix{i,:} = repmat(B{i},1,Q);
tempH{i,:} = IW{i,:} * Pn_train + BiasMatrix{i};
H{i,:} = 1 ./ (1 + exp(-tempH{i,:}));
A{i,:}=cond(H{i,:});
LW0{i,:} =pinv(H{i,:}') * Tn_train';
t{i,:}  = elmpredict(Pn_test,IW{i,:},B{i,:},LW0{i,:},'sig',0);
T{i,:} =mapminmax('reverse',t{i,:},ou);
L= length(P_test);
rmse1(i,:)=sqrt(sum((T_test-T{i,1}).^2)/L);

t2{i,:}  = elmpredict(Pn_train,IW{i,:},B{i,:},LW0{i,:},'sig',0);
T2{i,:} =mapminmax('reverse',t2{i,:},ou);
L= length(Pn_train);
rmse2(i,:)=sqrt(sum((T_train-T2{i,1}).^2)/L);
rmse(i,:)=mean(rmse1+rmse2);

r(i,:)=(L*sum(T{i,:}.*T_test)-sum(T{i,:})*sum(T_test))/sqrt(((L*sum((T{i,:}).^2)-(sum(T{i,:}))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
fitness=rmse2;
R2(i,:)=r(i,:);
pbest(i,:)=pop(i,:);
end

pFit = fitness;    
[ fMin, bestI ] = min( fitness );    % fMin denotes the global optimal adaptation value
bestX = pop( bestI, : );             % bestX denotes the global optimum
ffMin=fMin;

Lmin=1;
%% QPSO 
for a=1:q-1
    w=0.9-(0.9-0.4)*a/q;
    for i2=1:m
    %% Update particle velocity and position
    psi_i = exp(-(pop(i2,a) - xi).^2 / (2 * sigma^2)) / sqrt(2 * pi * sigma^2);
    psi_j = exp(-(bestX - xi).^2 / (2 * sigma^2)) / sqrt(2 * pi * sigma^2);
    delta_Psi = psi_j - psi_i;
    distance = norm(pop(i2,a) - bestX);
    pop(i2,a+1)=round(pop(i2,a)+w*pop(i2)+ beta * (delta_Psi / distance));
    nan_indices = isnan(pop);
    pop(i2,a+1) = round(40*rand(1));
    if pop(i2,a+1)<Lmin
        pop(i2,a+1)=Lmin
    end
        IW{i2,a+1} = rand(pop(i2,a+1),R) * 2 - 1;
        B{i2,a+1}= rand(pop(i2,a+1),1);
        BiasMatrix{i2,a+1} = repmat(B{i2,a+1},1,Q);
        tempH{i2,a+1} = IW{i2,a+1} * Pn_train + BiasMatrix{i2,a+1};
        H{i2,a+1} = 1 ./ (1 + exp(-tempH{i2,a+1}));
        A{i2,a+1}=cond(H{i2,a+1});
        LW0{i2,a+1} =pinv(H{i2,a+1}') * Tn_train';
        t{i2,a+1} = elmpredict(Pn_test,IW{i2,a+1},B{i2,a+1},LW0{i2,a+1},'sig',0);
        T{i2,a+1}=mapminmax('reverse',t{i2,a+1},ou);
        L= length(P_test);
        rmse1(i2,a+1)=sqrt(sum((T_test-T{i2,a+1}).^2)/L);
        r1(i2,a+1)=(L*sum(T{i2,a+1}.*T_test)-sum(T{i2,a+1})*sum(T_test))/sqrt(((L*sum((T{i2,a+1}).^2)-(sum(T{i2,a+1}))^2)*(L*sum((T_test).^2)-(sum(T_test))^2)));
       t2{i2,a+1}  = elmpredict(Pn_train,IW{i2,a+1},B{i2,a+1},LW0{i2,a+1},'sig',0);
T2{i2,a+1} =mapminmax('reverse',t2{i2,a+1},ou);
L= length(Pn_train);
rmse2(i2,a+1)=sqrt(sum((T_train-T2{i2,a+1}).^2)/L);
rmse(i2,a+1)=mean(rmse1(i2,a+1)+rmse2(i2,a+1)); 
        %% Individual location updates
        if rmse2(i2,a+1)<pFit(i2)
        pFit(i2)=rmse2(i2,a+1);
        pbest(i2)=pop(i2,a+1);
        end

        %% Collective position update
        if rmse2(i2,a+1)<ffMin 
        ffMin=rmse2(i2,a+1);
        bestX=pop(i2,a+1);
        end    
    end
        yyy(a)=ffMin;    
end  
    
plot(yyy)

[minValue, index] = min(rmse(:));
% Convert one-dimensional indexes to two-dimensional indexes
L= length(P_test);
[row, col] = ind2sub(size(rmse), index)
C=size(B{row,col},1)
t_F = elmpredict(Pn_test,IW{row,col},B{row,col},LW0{row,col},'sig',0);
T_F=mapminmax('reverse',t_F,ou);
rmse_F=sqrt(sum((T_test-T_F).^2)/L)
r_F=sqrt((L*sum(T_F.*T_test)-sum(T_F)*sum(T_test))/sqrt(((L*sum((T_F).^2)-(sum(T_F))^2)*(L*sum((T_test).^2)-(sum(T_test))^2))))
MAE_F=mean(abs(T_test - T_F))
MAPE_F=mean(abs((T_test - T_F) ./ T_test)) * 100    

save ('QPSO-ELM.MAT')