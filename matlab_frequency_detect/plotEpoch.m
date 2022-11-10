function [h,tle]=plotEpoch(EpochData_all,num,opt,fs,opt_seg,StEpochTime,chName,featShow,featN,flag_ncs,flag_acc)
 
EpochData=reshape(EpochData_all(num,:,:),size(EpochData_all,2),size(EpochData_all,3));
EpochData=detrend(EpochData);
EpochData=normalize(EpochData,1);

St=StEpochTime(num);

[feat_all,pkidx_all]=FeatInEpoch(EpochData,fs,opt);



h(1)=figure;
nfig=size(EpochData,2);
sz=9;

t=(1:length(EpochData))/fs;
cN={'green','red','red','red','blue','blue','blue'};
for i=1:nfig
    subplot(nfig,1,i);
    plot(t,EpochData(:,i),'LineWidth',0.5,'color',cN{i});
    hold on
    pkMax=pkidx_all{i}.max;pkMin=pkidx_all{i}.min;
    plot(t(pkMax),EpochData(pkMax,i),'^',...
         t(pkMin),EpochData(pkMin,i),'v');

    xlabel('time (s)','FontSize',sz)
    ylabel(chName{i},'FontSize',sz)
    
    
    txt_f=[];
    for j=1:length(featShow)  % for loop display all shown features 
      txt_f{j}=[' ',featN{featShow(j)},': ',num2str(feat_all(featShow(j),i),2)];
    
    end 
    
    title(join(txt_f),'FontSize',sz)

end

txt_flag=['NCS sel:',num2str(flag_ncs(num)),' ACC sel:',num2str(flag_acc(num))];
txt=['subj:', opt_seg.subjName,' seg:',opt_seg.segName,' side:',opt_seg.side,' Start T(min):',St/60];
tle=txt;
sgtitle([join(txt_flag) join(txt)],'fontsize', sz);
set(gcf,'Position',[100,10,800,900]);
set(gca,'fontsize', sz);





end