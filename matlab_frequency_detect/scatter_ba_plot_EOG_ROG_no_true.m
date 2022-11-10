close all


ExpDate='EOG_2';
i=4;
% ExpDate='EOG';
% i=1;
CaseName=['Case',num2str(i)];
fileName=[CaseName,'Routine1'];

LoadPath=['D:\eye RMG\data\','fig_REM_paper\','matFile\'];
SavePath=['D:\eye RMG\data\','fig_REM_paper\'];
load([LoadPath,ExpDate,fileName,'result_EMfreq.mat']);





%% plot scatter figure 
 a=linspace(0,70,20);
 b=a;
 sz=13;
 n=2;m=1;
 % only use 
ind=find(EM_f~=0);
EM_f=EM_f(ind);
br_eog=br_eog(ind);
br_rog=br_rog(ind);


Mean1 = mean(br_rog - br_eog);


[r1,p1] = corrcoef(br_rog,br_eog); 
r1 = r1(1,2);

h(1)=figure;
subplot(m,n,1)
% peak location 
scatter(br_rog ,br_eog,'x','MarkerEdgeColor','g');
hold on

plot(a,b,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2)
xlabel('ROG ','FontSize',sz)
ylabel('EOG','FontSize',sz)

max=60;
xlim([0 max]);
ylim([0 max]);

text(max*0.2,max*0.9,['\bf', 'r=',num2str(r1,2)],'FontSize',sz);




subplot(m,n,2)
Mean1 = mean(br_rog - br_eog);

std1 = std(br_rog - br_eog);

StdLim1 = [std1*1.96+Mean1, -std1*1.96+Mean1];


scatter((br_rog+br_eog)/2,(br_rog-br_eog)/2,'x','MarkerEdgeColor','g');
hold on


tt=30;
plot(a,Mean1.*ones(length(a),1),'color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth',1);
plot(a,StdLim1(1).*ones(length(a),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',1);
plot(a,StdLim1(2).*ones(length(a),1),'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',1);
text(tt,Mean1-0.2,['\bf', num2str(Mean1,2)],'FontSize',sz)
text(tt,StdLim1(1)-0.3,['\bf', num2str(StdLim1(1),2)],'FontSize',sz)
text(tt,StdLim1(2)+0.3,['\bf', num2str(StdLim1(2),2)],'FontSize',sz)




xlim([0 max]);
ylim([-6.5 6]);
xlabel('Mean ','FontSize',sz);
ylabel('Difference ','FontSize',sz);
sgtitle('EM rate (BPM)','FontSize',sz);
% title('wearable NCS  Forest','FontSize',sz);




set(gcf,'Position',[200,200,600,280]);

figName = [SavePath,ExpDate,fileName,'scatter_ba_eog_rog_no_true'];
print(h(1),[figName,'.tiff'],'-dtiff','-r300');
savefig(h(1),[figName,'.fig']);






