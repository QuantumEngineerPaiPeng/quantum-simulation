% 20180917
% Combine saved .fig files into a single one with multiple subplots
% https://www.mathworks.com/matlabcentral/answers/85219-subplot-of-already-saved-figures
% Load saved figures
% c1=[hgload('MedTauz4.fig');hgload('Medtauz3.fig');hgload('MedH3.fig');hgload('MedH5.fig')];
% c1=[hgload('Medtauz6.fig');hgload('Medtauz7.fig');hgload('MedH6.fig');hgload('MedH7.fig')];
c1=[hgload('MedH_110_4');hgload('MedH_110_5');hgload('Ztauz_lbits6_lbits_int');hgload('Ztauz_lbits6_lbits_int3')];
% c1=[hgload('dephase_pbitlbit.fig'),hgload('ll13.fig'),hgload('ll13.fig')];
% c1=[hgload('dephase_pbitlbit3.fig')];
% c1=[hgload('DisRealRare.fig'),hgload('MedZtauz_rare.fig')];
% c1=[hgload('MBLEM.fig'),hgload('MBLEM_BIC.fig'),hgload('MBLEM4cluster.fig')];
% c1=[hgload('Isingkmeans.fig'),hgload('IsingEM.fig')];
row=2;
col=2;
% Prepare subplots
figure
for p=1:length(c1)
    hh(p)=subplot(row,col,p);
    axtemp=get(c1(p),'CurrentAxes');
    copyobj(allchild(axtemp),hh(p));
    hh(p).XScale=axtemp.XScale;
    hh(p).YScale=axtemp.YScale;
    hh(p).XLim=axtemp.XLim;
    hh(p).YLim=axtemp.YLim;
    hh(p).XLabel=axtemp.XLabel;
    hh(p).YLabel=axtemp.YLabel;
    ppStyle(20,2,10)
end
