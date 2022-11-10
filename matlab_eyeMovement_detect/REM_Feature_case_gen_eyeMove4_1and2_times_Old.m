%data normalized to N(0,1)
%devide into epochs  


close all 
clear all
% waveform left shift 'a': time +a; waveform right shift 'a': time -a; 


opt.delay=zeros(34,1);


% ind=1:16; 
% opt.delay(ind)=opt.delay(ind)+0.5;

opt.delay=opt.delay-1;

caseNum=2;

filt=[0.1,5];
for i=caseNum
  
  main (i,opt,filt)
end
filt=[0.05,10];
for i=caseNum
  
  main (i,opt,filt)
end

function main(i,opt,filt)
close all
ExpDate='5_10';
dataPath=['D:\eye RMG\data\',ExpDate,'\'];
SavePath=['D:\eye RMG\data\',ExpDate,'\','feature\'];

Feat_ver=['filt_',num2str(filt(1)),'_',num2str(filt(2))];



ch_plot1=1;ch_plot2=31;
CaseName=['Case',num2str(i)];
opt.constantPad=1;

fileName=[CaseName,'Routine4'];
fs=5e3;
fsDS=500;
toff=[5:170]';



period=[toff(1),toff(end)];

% if exist([SavePath,'\',Feat_ver,'\',CaseName,'.mat'],'file')
%     fprintf('found  case %s\n',fileName);
%     return
% end

filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];

if ~exist(filePathName_m,'file')
   convertTDMS(true,filePathName);
end
load(filePathName_m);
gesture={,'Up','Up','Up','Down','Down','Down','Down',...
     'Right','Right','Right','Right','Left','Left','Left','Left','0',...
    'Up×2','Up×2','Up×2','Up×2','Down×2','Down×2','Down×2','Down×2',...
    'Right×2','Right×2','Right×2','Right×2', 'Left×2','Left×2', 'Left×2','Left×2'...
    };
gesture_all={''};
for i=1:length(gesture)
gesture_all{i+1}=gesture{i};
end    
    
load([dataPath,fileName,'.mat']);
Chan_Name={'Tx1Rx1 amp','TxRx1 ph','Tx2Rx1 amp','Tx2Rx1 ph','Tx3Rx1 amp','Tx3Rx1 ph','Tx4Rx1 amp','Tx4Rx1 ph',...
    'Tx1Rx2 amp','Tx1Rx2 ph','Tx2Rx2 amp','Tx2Rx2 ph','Tx3Rx2 amp','Tx3Rx2 ph','Tx4Rx2 amp','Tx4Rx2 ph',...
    'Tx1Rx3 amp','TxRx3 ph','Tx2Rx3 amp','Tx2Rx3 ph','Tx3Rx3 amp','Tx3Rx3 ph','Tx4Rx3 amp','Tx4Rx3 ph',...
    'Tx1Rx4 amp','TxRx4 ph','Tx2Rx4 amp','Tx2Rx4 ph','Tx3Rx4 amp','Tx3Rx4 ph','Tx4Rx4 amp','Tx4Rx4 ph'}';

Ch_num=cat(1,[3:18]',[23:38]');


% filter type 
opt.filtType = 'LpHp'; opt.orderHP = 5;
opt.f3db = filt(1); opt.fpLP = filt(2); opt.fstLP = opt.fpLP+1;

for i = 1:32
  
  Ch_data(:,i)=ConvertedData.Data.MeasuredData(Ch_num(i)).Data; % ch1 amp
  Ch_data_raw=Ch_data;
    if mod(i,2) == 0  % if phase channel 
      Ch_data(:,i)=unwrap(deg2rad(Ch_data(:,i)));
  end
%   
  [Ch_data_filt(:,i),Ch_data_filt_raw(:,i)]= filtSignal(Ch_data(:,i),toff,opt);
  Ch_data_filt(:,i)=detrend(Ch_data_filt(:,i));
end 


% 16 channels do vector normalization to get real and imagination part
opt.Vec_filtType='lowpass';opt.Vec_fLow=0.05;opt.Vec_fHigh=5;
for j=1:16
    
   [Ch_data_complex(:,j)]=vector_norm(Ch_data_raw(:,2*j-1),Ch_data_raw(:,2*j),toff,opt,fsDS,fs);


end

%% waveform segment into epochs 
% start from 5s,

 
label=[1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,0,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8];

tWin=5; tSlide=5; % no padding 
tStart=5;

% start / end time of each epoch  2* epochNum
tEnd=toff(end)-toff(1)-tStart; 
StartEndT=cat(1,tStart:tSlide:tEnd, tStart+tWin:tSlide:tEnd+tWin);  


%  special time stamp add time delay
for i=1:length(StartEndT)


StartEndT(:,i)=StartEndT(:,i)+opt.delay(i); 

end






% divide waveform into epochs 
% NOT PADDING 5S SEGMENT  directly 
for i= 1:length(label)

Ch_data_epoch(:,:,i)=Ch_data_filt((StartEndT(1,i)*fsDS)+1:StartEndT(2,i)*fsDS,:);

Ch_data_epoch_complex(:,:,i)=Ch_data_complex((StartEndT(1,i)*fsDS)+1:StartEndT(2,i)*fsDS,:);



end

%% plot figures 


%plot whole waveforms 4 main channels 
h_w=plotWaveform(Ch_data_filt,CaseName,gesture_all,fsDS);

h_e1=plotEpoch(Ch_data_epoch,ch_plot1,gesture,fsDS,StartEndT,Chan_Name,fileName);
h_e2=plotEpoch(Ch_data_epoch,ch_plot2,gesture,fsDS,StartEndT,Chan_Name,fileName);



for i=1:length(h_w)
figName = [dataPath,'fig_case','\',fileName,'_waveform_',num2str(i)];
print(h_w(i),[figName,'.tiff'],'-dtiff','-r300');
savefig(h_w(i),[figName,'.fig']);
end
for j=1:length(h_e1)
figName = [dataPath,'fig_case','\',fileName,'ch',num2str(ch_plot1),'_Epoch_',num2str(j)];
print(h_e1(j),[figName,'.tiff'],'-dtiff','-r300');

end
for j=1:length(h_e2)
figName = [dataPath,'fig_case','\',fileName,'ch',num2str(ch_plot2),'_Epoch_',num2str(j)];
print(h_e2(j),[figName,'.tiff'],'-dtiff','-r300');

end

% delete some 0 epochs 

delind=[16]; % delete 
Ch_data_epoch(:,:,delind)=[];
Ch_data_epoch_complex(:,:,delind)=[];
label_all=label; % all labels without delete 
label(delind)=[];
gesture(delind)=[];

opt.delind=delind;opt.StartEndTime=StartEndT;opt.gesture_all=gesture_all;opt.label_all=label_all;
opt.dataPath=filePathName_m;opt.Case=fileName;opt.toff=toff;

save([SavePath,'\',Feat_ver,'\',fileName,'.mat'],'Ch_data_epoch','Ch_data_epoch_complex','label','gesture','Ch_data_complex','Ch_data_filt','Chan_Name','opt');
end



function [ampCh_filt_norm,ampCh_filt]=filtSignal(ampCh,toff,opt)
    
    fs=5e3;
    fsDS=500;
    
    

    ampCh=resample(ampCh,fsDS,fs);
    

    ampCh=ampCh((toff(1)*fsDS):toff(size(toff))*fsDS);
    
  

    ampCh_filt = filterLpHp(ampCh,fsDS,opt); % th amp
    
%     [,PS] = mapminmax(ampCh_filt');
%     ampCh_filt_norm=ampCh_filt_norm';
    ampCh_filt_norm = normalize(ampCh_filt);

end 


