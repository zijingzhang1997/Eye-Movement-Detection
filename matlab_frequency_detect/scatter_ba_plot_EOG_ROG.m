close all


ExpDate='EOG_2';
i=4;
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


Mean1 = mean(EM_f - br_eog);
Mean2 = mean(EM_f - br_rog);

[r1,p1] = corrcoef(EM_f,br_eog); 
r1 = r1(1,2);
[r2,p1] = corrcoef(EM_f,br_rog); 
r2 = r2(1,2);

h(1)=figure;
subplot(m,n,1)
% peak location 
scatter(EM_f ,br_eog,'x','MarkerEdgeColor','b');
hold on
scatter(EM_f ,br_rog,'*','MarkerEdgeColor','r');
plot(a,b,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2)
xlabel('True ','FontSize',sz)
ylabel('Estimation','FontSize',sz)

max=60;
xlim([0 max]);
ylim([0 max]);
legend('EOG','ROG',...
'FontSize',sz,'location','northwest','orientation','horizontal');
legend('boxoff')
text(max*0.2,max*0.9,['\bf', 'r_{ROG}=',num2str(r2,2)],'FontSize',sz);
text(max*0.2,max*0.8,['\bf', 'r_{EOG}=',num2str(r1,2)],'FontSize',sz);



subplot(m,n,2)
Mean1 = mean(EM_f - br_eog);
Mean2 = mean(EM_f - br_rog);
std1 = std(EM_f - br_eog);
std2 = std(EM_f - br_rog);
StdLim1 = [std1*1.96+Mean1, -std1*1.96+Mean1];
StdLim2 = [std2*1.96+Mean2, -std2*1.96+Mean2]; 

scatter((EM_f+br_eog)/2,(EM_f-br_eog)/2,'x','MarkerEdgeColor','b');
hold on
scatter((EM_f+br_rog)/2,(EM_f-br_rog)/2,'*','MarkerEdgeColor','r');

tt=30;
plot(a,Mean1.*ones(length(a),1),'color','b','LineStyle','-.','LineWidth',1);
plot(a,StdLim1(1).*ones(length(a),1),'color','b','LineStyle',':','LineWidth',1);
plot(a,StdLim1(2).*ones(length(a),1),'color','b','LineStyle',':','LineWidth',1);
text(tt,Mean1-0.2,['\bf', num2str(Mean1,2)],'FontSize',sz,'color','b')
text(tt,StdLim1(1)-0.3,['\bf', num2str(StdLim1(1),2)],'FontSize',sz,'color','b')
text(tt,StdLim1(2)+0.3,['\bf', num2str(StdLim1(2),2)],'FontSize',sz,'color','b')


plot(a,Mean2.*ones(length(a),1),'color','r','LineStyle','-.','LineWidth',1);
plot(a,StdLim2(1).*ones(length(a),1),'color','r','LineStyle',':','LineWidth',1);
plot(a,StdLim2(2).*ones(length(a),1),'color','r','LineStyle',':','LineWidth',1);
text(tt,Mean2-0.2,['\bf', num2str(Mean2,2)],'FontSize',sz,'color','r')
text(tt,StdLim2(1)-0.3,['\bf', num2str(StdLim2(1),2)],'FontSize',sz,'color','r')
text(tt,StdLim2(2)+0.3,['\bf', num2str(StdLim2(2),2)],'FontSize',sz,'color','r')

xlim([0 max]);
xlabel('Mean ','FontSize',sz);
ylabel('Difference ','FontSize',sz);
sgtitle('EM rate (BPM)','FontSize',sz);
% title('wearable NCS  Forest','FontSize',sz);




set(gcf,'Position',[200,200,600,280]);

figName = [SavePath,ExpDate,fileName,'scatter_ba_eog_rog'];
print(h(1),[figName,'.tiff'],'-dtiff','-r300');
savefig(h(1),[figName,'.fig']);






