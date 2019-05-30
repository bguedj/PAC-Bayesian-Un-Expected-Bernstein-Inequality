close all
%% Loading the data 
prefix = "_sigmoid-synthetic-imformedPrior-FALSE-";
dim = 10;
name = prefix+num2str(dim)+".mat";

load("save/results"+name)
results = labpcexport;
BProb = mean(results(:,1,:),1); BProb = BProb(:);
BTS = mean(results(:,2,:),1); BTS = BTS(:);
BKL = mean(results(:,3,:),1); BKL = BKL(:);
BCT = mean(results(:,4,:),1); BCT = BCT(:);

load("save/Vn"+name)
Vn = labpcexport;
MVn = mean(Vn,1);

load("save/VnPrim"+name)
Vnprim = labpcexport;
MVnprim = mean(Vnprim,1);

load("save/LnERMtrain"+name)
LnERMtrain = labpcexport;
MLnERM = mean(LnERMtrain,1);


%% Constants
lw = 1.5;
f=1;

left= 1;
bottom =1;
width =3;
height=6;

%% Plotting the left plot
F1 = figure('units','inches','position',[0 0 2*width+2 height+1.6]);

i=1;
l = (i-1) * width + left;
b = bottom;

h(i) = axes('units','inches','Position',[l b width height]);
hold on 
box on
plot(1:10, BProb, 'b*-','linewidth',lw)
plot(1:10, BKL, 'r^-','linewidth',lw)
plot(1:10, BTS, 'go-','linewidth',lw)
plot(1:10, BCT, 'ks-','linewidth',lw)
plot(1:10, MVn, 'k:','linewidth',lw*f)
plot(1:10, MVnprim, 'k--','linewidth',lw*f)
plot(1:10, MLnERM, 'k-','linewidth',lw*f)


axis([0.5 10.5 0.0 0.5])
set(h(i), 'fontsize', 20,'XTick', linspace(1.25,8.75,4), 'XTickLabel', {1000,3000,5000,7000,}) 
xlabel('SAMPLE SIZE')
ylabel('BOUND')
title("Dimension = " + num2str(dim))
lgd =legend('Our Bound', 'Maurer','TS','Catoni',"$V_n$","$V'_n$","$L_n(\hat{h}(Z_{\leq n}))$");
set(lgd, 'Interpreter','latex','box','off','FontSize',22)

set(gca,'FontName', 'Times New Roman')

%% Plotting the right plot
dim=50;
name = prefix+num2str(dim)+".mat";

load("save/results"+name)
results = labpcexport;
BProb = mean(results(:,1,:),1); BProb = BProb(:);
BTS = mean(results(:,2,:),1); BTS = BTS(:);
BKL = mean(results(:,3,:),1); BKL = BKL(:);
BCT = mean(results(:,4,:),1); BCT = BCT(:);

load("save/Vn"+name)
Vn = labpcexport;
MVn = mean(Vn,1);

load("save/VnPrim"+name)
Vnprim = labpcexport;
MVnprim = mean(Vnprim,1);

load("save/LnERMtrain"+name)
LnERMtrain = labpcexport;
MLnERM = mean(LnERMtrain,1);


i=2;
l = (i-1) * width + left;
b = bottom;
h(i) = axes('units','inches','Position',[l b width height]);
hold on 
box on
plot(1:10, BProb, 'b*-','linewidth',lw)
plot(1:10, BKL, 'r^-','linewidth',lw)
plot(1:10, BTS, 'go-','linewidth',lw)
plot(1:10, BCT, 'ks-','linewidth',lw)
plot(1:10, MVn, 'k:','linewidth',lw*f)
plot(1:10, MVnprim, 'k--','linewidth',lw*f)
plot(1:10, MLnERM, 'k-','linewidth',lw*f)


axis([0.5 10.5 0.0 0.5])
set(h(i), 'fontsize', 20,'XTick', linspace(1.25,8.75,4), 'XTickLabel', {1000,3000,5000,7000,}) 
xlabel('SAMPLE SIZE')
title("DIMENSION = " + num2str(dim))
%lgd =legend('Our Bound', 'Maurer','Tolstikhin','Catoni',"$V_n$","$V'_n$","$L_n(\hat{h}(Z_{\leq n}))$");
%set(lgd, 'Interpreter','latex','box','off', 'FontSize',22)
set(h(i),'yticklabel',[])
set(gca,'FontName', 'Times New Roman')

%% Saving the Figure
saveas(gcf,'Figures/synthetic_exp2','epsc')
%close all


