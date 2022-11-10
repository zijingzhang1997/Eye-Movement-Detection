
clear all
close all
i=1;
dataPath=['D:\eye RMG\data\4_9\'];
SavePath=[dataPath,'fig_case'];
CaseName=['Case',num2str(i)];
fileName=[CaseName,'Routine1'];
fs=5e3;
fsDS=500;
toff=[5:58]';
period=[toff(1),toff(end)];
msize=30;


filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];
if ~exist(filePathName_m,'file')
   convertTDMS(true,filePathName);
end
load(filePathName_m);


load([dataPath,fileName,'.mat']);

opt.filtType = 'LpHp'; opt.orderHP = 5;
opt.f3db = 0.05; opt.fpLP = 20; opt.fstLP = opt.fpLP+1;
Ch_num=cat(1,[3:18]',[23:38]');
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

%%
opt_emg.filtType = 'LpHp'; opt_emg.orderHP = 5;
opt_emg.f3db = 0.5; opt_emg.fpLP = 15; opt_emg.fstLP =opt_emg.fpLP+1;
opt_emg.SMorder=3;opt_emg.SMlen=99;  %% smoothing filter parameter , order less, length longer, more smooth 
opt_emg.movmean=200;
% opt_emg.SMlen=199; opt_emg.movmean=100;
Ch_num_emg=[19,20];
for i = 1:2  %% use first two channels of EMG 
  
  Ch_data_emg(:,i)=ConvertedData.Data.MeasuredData(Ch_num_emg(i)).Data; % ch1 amp
  
  [Ch_data_emg_filt(:,i),Ch_data_emg_filt_raw(:,i),Ch_data_emg_en(:,i),Ch_data_emg_raw(:,i)]= filtSignal_emg(Ch_data_emg(:,i),toff,opt_emg);
  % Ch_data_emg_filt_raw : only filter, no smooth 
  % Ch_data_emg_filt  : filter + smooth + normalization
  Ch_data_emg_filt(:,i)=detrend(Ch_data_emg_filt(:,i));
end 
opt_emg.smOpt=0;
if opt_emg.smOpt==1
    Ch_data_emg_filt=smoothdata(Ch_data_emg_filt,1,'movmean',opt_emg.movmean);
end

Ch_data_filt=[Ch_data_filt,Ch_data_emg_filt];
Ch_data_filt=detrend(Ch_data_filt);
Ch_data_emg=cat(3,Ch_data_emg_raw,Ch_data_emg_en,Ch_data_emg_filt_raw,Ch_data_emg_filt);




Chan_Name={'Tx1Rx1 amp','TxRx1 ph','Tx2Rx1 amp','Tx2Rx1 ph','Tx3Rx1 amp','Tx3Rx1 ph','Tx4Rx1 amp','Tx4Rx1 ph',...
    'Tx1Rx2 amp','Tx1Rx2 ph','Tx2Rx2 amp','Tx2Rx2 ph','Tx3Rx2 amp','Tx3Rx2 ph','Tx4Rx2 amp','Tx4Rx2 ph',...
    'Tx1Rx3 amp','TxRx3 ph','Tx2Rx3 amp','Tx2Rx3 ph','Tx3Rx3 amp','Tx3Rx3 ph','Tx4Rx3 amp','Tx4Rx3 ph',...
    'Tx1Rx4 amp','TxRx4 ph','Tx2Rx4 amp','Tx2Rx4 ph','Tx3Rx4 amp','Tx3Rx4 ph','Tx4Rx4 amp','Tx4Rx4 ph',...
    'EMG 1','EMG 2',}';


opt.pmin=0.3;opt.hmin=0.5;opt.dmin=0.05;opt.dmax=0.8;opt.dis=0.1;


ChPlot=[1,11,21,31];
h_all(1)=plotWaveform_nopeak(Ch_data_filt(:,ChPlot),Chan_Name(ChPlot),fsDS,msize);
ChPlot=[2,12,22,32];
h_all(2)=plotWaveform_nopeak(Ch_data_filt(:,ChPlot),Chan_Name(ChPlot),fsDS,msize);

% SavePath='D:\RFMG\matlab code\for_fig_manuscript\fig\timeTest\';

figName = [SavePath,'\',fileName,'Waveform_amp'];
print(h_all(1),[figName,'.tiff'],'-dtiff','-r300');
savefig(h_all(1),[figName,'.fig']);
figName = [SavePath,'\',fileName,'Waveform_ph'];
print(h_all(2),[figName,'.tiff'],'-dtiff','-r300');
savefig(h_all(2),[figName,'.fig']);

%  for i=3 can detect all peaks opt.p1min=0.2;opt.h1min=0.5;  RMG 31
%  opt.p2min=0.7;opt.h2min=0.3; EMG 34

% for i=4 can detect all peaks opt.p1min=0.2;opt.h1min=0.2; opt.p2min=0.7;opt.h2 min=0.3;




[p1,loc1,wid1,pro1]=findpeaks(Ch_data_filt(:,11),fsDS,'MinPeakProminence',opt.pmin,'MinPeakHeight',opt.hmin,...
    'MinPeakWidth',opt.dmin,'MaxPeakWidth',opt.dmax);


[p2,loc2,wid2,pro2]=findpeaks(Ch_data_filt(:,34),fsDS,...
    'MinPeakProminence',opt.pmin,'MinPeakHeight',opt.hmin,'MinPeakWidth',opt.dmin,'MaxPeakWidth',opt.dmax);

interT_1=diff(loc1);
interT_2=diff(loc2);
inT1=[mean(interT_1) std(interT_1)];
inT2=[mean(interT_2) std(interT_2)];

pro1m=[mean(pro1) std(pro1)];
pro2m=[mean(pro2) std(pro2)];

%%
function h=plotWaveform_s(Ch_data,ChName,fsDS,msize,opt)

tOff=((0:(length(Ch_data)-1))/fsDS)';

h(1)=figure;
a=0;
sz=11;
cN={'blue','red'};
n=size(Ch_data,2);
for i =1:n
subplot(n,1,i);

plot(tOff,Ch_data(:,i),'LineWidth',0.5,'color',cN{i},'LineWidth',1);
[pk,loc,wid2,pro2]=findpeaks(Ch_data(:,i),fsDS,...
    'MinPeakProminence',opt.pmin,'MinPeakHeight',opt.hmin,...
    'MinPeakWidth',opt.dmin,'MaxPeakWidth',opt.dmax);
hold on

scatter(tOff(round(loc*fsDS)),pk,msize,'green','filled','^')
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([a round(max(tOff))])

legend(ChName{i},'FontSize',sz,'Location','northoutside')
legend('boxoff')
end

set(gcf,'Position',[100,100,450,300]);






end 
function h=plotWaveform_nopeak(Ch_data,ChName,fsDS,msize)

tOff=((0:(length(Ch_data)-1))/fsDS)';

h(1)=figure;
a=0;
sz=13;
cN={'blue','red','green','m'};
n=size(Ch_data,2);

for i =1:n
subplot(n,1,i);

plot(tOff,Ch_data(:,i),'LineWidth',0.5,'color',cN{i},'LineWidth',1);

xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([a round(max(tOff))])

legend(ChName{i},'FontSize',sz,'Location','northoutside')
legend('boxoff')
end

set(gcf,'Position',[100,100,1500,1000]);






end 

function [ampCh_filt_norm,ampCh_filt]=filtSignal(ampCh,toff,opt)
    
    fs=5e3;
    fsDS=500;
    
  
    ampCh=resample(ampCh,fsDS,fs);
    ampCh=ampCh((toff(1)*fsDS):toff(size(toff))*fsDS);
    ampCh_filt = filterLpHp(ampCh,fsDS,opt); % th amp
    
    ampCh_filt_norm = normalize(ampCh_filt);

end 

function [ampCh_filt_norm,ampCh_filt,ampCh_en,ampCh_raw]=filtSignal_emg(ampCh,toff,opt)
    
    fs=1e3;  % biopac sampling rate =1e3 
    fsDS=500; % down sampling rate  same as NCS  
    
  
    ampCh=resample(ampCh,fsDS,fs);
    ampCh_raw=ampCh((toff(1)*fsDS):toff(size(toff))*fsDS);
    
    [ampCh_en,~] = envelope(ampCh_raw,50,'peak');
    
    ampCh_en=detrend(ampCh_en);
    ampCh_filt = filterLpHp(ampCh_en,fsDS,opt); % th amp
    
    
    ampCh_filt_sm = sgolayfilt(ampCh_filt,opt.SMorder,opt.SMlen);
    %ampCh_filt_sm = smoothdata(ampCh_filt,'movmean',200);
    
    ampCh_filt_norm = normalize(ampCh_filt_sm);
    
end 

