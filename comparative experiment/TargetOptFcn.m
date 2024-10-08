function targetValue = TargetOptFcn(r,U,C,D)
% Input:  Name     description         type
%          r         属性编码          array
%          U           论域            array
%          C         条件属性          array
%          D         决策属性          array
% Output: Name      description        tpye
%      targetValue    目标值           double

%==========================================================================

    if(length(r)~= size(C,2))
        error('RS:error','wrong binary input r.');
    end
    r_bin =round(r);
    C_reduct = C(:,r_bin==1);
    % 第一个目标是最大化目标D对于C_reduct的依赖度
    target1 = DependencyDegree(U,C_reduct,D);
    % 第二个目标是尽量减小属性C中1的个数，也就是尽量的约简C
    target2 = (length(r_bin) - sum(r_bin))/length(r_bin);
    % 第三个目标是让染色体Lr中的数字尽可能的接近1或者0
    target3 = sum(abs(r.*(r-1)))/length(r);
    % 注意这里为了更好的优化结果多3个目标进行了加权处理，实际上target1 +target2
    % 就是论文里的目标函数
    % ga函数默认是最小化目标函数，所以这里前2个目标要加个负号转为最小化
    targetValue = -5*target1 - 2*target2 + target3;   
end
