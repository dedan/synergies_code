function figure1 = createfigure(YMatrix1)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 01-Dec-2008 15:38:01

% Create figure
figure1 = figure('XVisual',...
   '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)','Visible','off');

% Create axes
axes1 = axes('Parent',figure1);
hold('all');

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1);
set(plot1(1),'DisplayName','nmf pronation');
set(plot1(2),'DisplayName','nmf suplination');
set(plot1(3),'DisplayName','90 %','LineStyle','--');
set(plot1(4),'DisplayName','pcaica pronation');
set(plot1(5),'DisplayName','pcaica suplination');


% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.6437 0.5265 0.174 0.1439]);


