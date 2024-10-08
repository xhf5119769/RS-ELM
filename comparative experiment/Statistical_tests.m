clc
clear

load('RS-ELM.mat')      
for i=1:6
          LWi=V(:,i)*s(i)*U(:,i)'*Tn_train';
          LW_WT=LWi+LW_WT;
          BiasMatrix_WT = repmat(B2,1,22);
          tempH_WT = IW2 * Pn_test + BiasMatrix_WT;
          H_WT2 =W*( 1 ./ (1 + exp(-tempH_WT)));
          t_WT = (H_WT2' * LW_WT)';
          T_WT = mapminmax('reverse',t_WT,ou);
E1=(T_test-T_WT);
      end
data1=E1';

load('IR-ELM.mat')
IW=IWW(1:28,:);
B=BB(1:28,:);
LW=LWW(1:28,:);
t_new = elmpredict(Pn_test,IW,B,LW,'sig',0);
T_new=mapminmax('reverse',t_new,ou);
E2=(T_test-T_new);
data2=E2';

load('GP-ELM.mat')
E3=(T_test-pruned_T);
data3=E3';

load('QPSO-ELM.mat')
E4=(T_test-T_F);
data4=E4';

data=[data1 data2 data3 data4];


[p, tbl, stats] = friedman(data, 1)

disp(['The p-value for the Friedman test is:', num2str(p)]);

% Determine whether it is significant
alpha = 0.05; % significance level
if p < alpha
    disp('Differences between models are statistically significant.');
else
    disp('Differences between models are not statistically significant.');
end

[c, m, h, gnames] = multcompare(stats, 'CType', 'dunn-sidak');

% Interpreting the results of multiple comparisons
for i = 1:size(c,1)
    group1 = gnames(c(i,1));
    group2 = gnames(c(i,2));
    lower_limit = c(i,3);
    mean_diff = c(i,4);
    upper_limit = c(i,5);
    p_value = c(i,6);
     % fprintf('模型 %s 与模型 %s 的比较：\n', group1, group2);
    % fprintf('  均值差异：%.4f\n', mean_diff);
    % fprintf('  95%%置信区间：[%.4f, %.4f]\n', lower_limit, upper_limit);
    % fprintf('  p值：%.4f\n', p_value);
        if p_value < alpha
        disp('  The difference is significant。');
    else
        disp('  The difference is not significant。');
    end
    disp('-----------------------');
end

