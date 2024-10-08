clc
clear
load('y.mat')
load('X_dosc.mat')

X=X_DOSC;
y=y;
%% original spectra
R=corr(X', y,'type','Pearson');
%% Reciprocal spectra
X_R=1./X;
R_R=corr(X_R', y,'type','Pearson');
%% Square root spectra
X_RMS=sqrt(X);
R_RMS=corr(X_RMS', y,'type','Pearson');
%% Logarithmic spectra
X_L=log(X);
R_L=corr(X_L', y,'type','Pearson');
%% First-order derivative spectra
% X_FD=diff(X);
% R_FD=corr(X_FD', y,'type','Pearson');

window_size = 3;  % 窗口大小（奇数）
order = 1;  % 阶数
sg_coeff = sgolay(order, window_size);
X_FD=[];
for i=1:82
RFD = conv(X(:,i), sg_coeff(:, order+1), 'same');
X_FD(:,i)=RFD;
end
R_FD=corr(X_FD', y,'type','Pearson');
%% Second-order derivative spectra
% X_SD=diff(X,2);
window_size = 3;  % 窗口大小（奇数）
order = 2;  % 阶数
sg_coeff = sgolay(order, window_size);
X_SD=[];
for i=1:82
RSD = conv(X(:,i), sg_coeff(:, order+1), 'same');
X_SD(:,i)=RSD;
end
R_SD=corr(X_SD', y,'type','Pearson');
%% OS
figure (1)
axes1 = axes('Parent',figure (1),...
    'Position',[0.139839034205231 0.143939393939394 0.823943661971832 0.83030303030303]);
plot(R, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OS');
hold(axes1,'on');
 box on; %开启右面和上面的坐标轴
ylabel('Pearson correlation coefficient');
xlabel('Band number');
set(axes1,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'B1','B2','B3','B4','B5','B6','B7','B8','B8a','B9','B11','B12'});
%% Reciprocal spectra
figure (2)
axes1 = axes('Parent',figure (2),...
    'Position',[0.139839034205231 0.143939393939394 0.823943661971832 0.83030303030303]);
hold(axes1,'on');
plot(R_R, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'RecS');
 box on; %开启右面和上面的坐标轴
ylabel('Pearson correlation coefficient');
xlabel('Band number');
set(axes1,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'B1','B2','B3','B4','B5','B6','B7','B8','B8a','B9','B11','B12'});
%% Square root spectra
figure (3)
axes1 = axes('Parent',figure (3),...
    'Position',[0.139839034205231 0.143939393939394 0.823943661971832 0.83030303030303]);
plot(R_RMS, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'SRS');
hold(axes1,'on');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
set(axes1,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'B1','B2','B3','B4','B5','B6','B7','B8','B8a','B9','B11','B12'});
%% Logarithmic spectra
figure (4)
axes1 = axes('Parent',figure (4),...
    'Position',[0.139839034205231 0.143939393939394 0.823943661971832 0.83030303030303]);
plot(R_L, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'LS');
hold(axes1,'on');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
set(axes1,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'B1','B2','B3','B4','B5','B6','B7','B8','B8a','B9','B11','B12'});
%% First-order derivative spectra
figure (6)
axes1 = axes('Parent',figure (6),...
    'Position',[0.139839034205231 0.143939393939394 0.823943661971832 0.83030303030303]);
hold(axes1,'on');
plot(R_FD, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'FODS');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
set(axes1,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'B1','B2','B3','B4','B5','B6','B7','B8','B8a','B9','B11','B12'});
%% Second-order derivative spectra
figure (7)
axes1 = axes('Parent',figure (7),...
    'Position',[0.139839034205231 0.143939393939394 0.823943661971832 0.83030303030303]);
hold(axes1,'on');
plot(R_SD, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'SODS');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
ylabel('Pearson correlation coefficient');
xlabel('Band number');
set(axes1,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'B1','B2','B3','B4','B5','B6','B7','B8','B8a','B9','B11','B12'});
%% Saving Variables
save('X_FD','X_FD')
