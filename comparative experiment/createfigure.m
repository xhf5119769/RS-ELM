function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  x 数据的向量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 15-May-2023 21:03:50 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1,'Position',[0.092 0.177 0.875 0.8]);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot(X1,YMatrix1,'Parent',axes1);

% 创建 ylabel
ylabel('Reflectance (%)');

% 创建 xlabel
xlabel(['Wavelength (nm)';'(d)            ']);

% 取消以下行的注释以保留坐标区的 X 范围
 xlim(axes1,[300 2600]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes1,[-Inf Inf]);
box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'FontSize',24,'XTick',[300 600 900 1200 1500 1800 2100 2400 2600]);
