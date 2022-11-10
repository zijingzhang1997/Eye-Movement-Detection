
clear all
close all
i=3;
ExpDate='EOG_2';
dataPath=['D:\eye RMG\data\',ExpDate,'\'];
addpath('D:\eye RMG\matlab code_all');
addpath('D:\COVID\COVID_HF_spectrum\dyspnea_study_new_code');

CaseName=['Case',num2str(i)];
fileName=[CaseName,'Routine1'];
fs=5e3;
fsDS=500;

filt=[0.1,2];
% make sure filter LF is lower so higher frequency is filtered out

% Less noisy for peak detection

%toff=[10:length(Ch_data_raw)/fs-2]';

filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];
if ~exist(filePathName_m,'file')
   convertTDMS(true,filePathName);
end
load(filePathName_m);


load([dataPath,fileName,'.mat']);

opt.filtType = 'LpHp'; opt.orderHP = 4;
opt.f3db = filt(1); opt.fpLP = filt(2); opt.fstLP = opt.fpLP+1;
Ch_num=cat(1,[3:18]',[23:38]');

% ground truth of movement 


t=1:170;
MotionTrue=zeros(length(t),1);
MotionTrue(1:28)=10;MotionTrue(33:63)=15;MotionTrue(33:63)=15;
MotionTrue(68:98)=20;MotionTrue(103:131)=30;MotionTrue(136:166)=60;
Motion=[t'  MotionTrue];

for i = 1:32
  
  Ch_data(:,i)=ConvertedData.Data.MeasuredData(Ch_num(i)).Data; % ch1 amp
  
  Ch_data_raw=Ch_data;
    if mod(i,2) == 0  % if phase channel 
      Ch_data(:,i)=unwrap(deg2rad(Ch_data(:,i)));
    end
  toff=[10:length(Ch_data_raw)/fs-2]';  
  [Ch_data_filt(:,i),Ch_data_filt_raw(:,i)]= filtSignal(Ch_data(:,i),toff,opt,fs,fsDS);
  Ch_data_filt(:,i)=detrend(Ch_data_filt(:,i));
end 

%%  EOG 
% EOG is less noisy, can directly use the same filter of low pass 
% so noisy peaks can be smoothed , better for peak detection 
opt_emg.filtType = 'LpHp'; opt_emg.orderHP = 5;
opt_emg.f3db = 0.1; opt_emg.fpLP = 2; opt_emg.fstLP =opt_emg.fpLP+1;
opt_emg.SMorder=5;opt_emg.SMlen=9;  %% smoothing filter parameter , order less, length longer, more smooth 
opt_emg.movmean=50;
opt_emg.fsDS=fsDS;
% opt_emg.SMlen=199; opt_emg.movmean=100;
Ch_num_emg=[19];
for i = 1  %% use first two channels of EMG 
  
  Ch_data_emg(:,i)=ConvertedData.Data.MeasuredData(Ch_num_emg(i)).Data; % ch1 amp
  
  [Ch_data_emg_filt(:,i),Ch_data_emg_filt_raw(:,i),Ch_data_emg_en(:,i),Ch_data_emg_raw(:,i)]= filtSignal_emg(Ch_data_emg(:,i),toff,opt_emg);
  % Ch_data_emg_filt_raw : only filter, no smooth 
  % Ch_data_emg_filt  : filter + smooth + normalization
  Ch_data_emg_filt(:,i)=detrend(Ch_data_emg_filt(:,i));
end 


Ch_data_filt=[Ch_data_filt,Ch_data_emg_filt];
Ch_data_filt = sgolayfilt(Ch_data_filt,4,51);
Ch_data_emg=cat(3,Ch_data_emg_raw,Ch_data_emg_en,Ch_data_emg_filt_raw,Ch_data_emg_filt);


%%  peak detection + feature extraction 
Ch_data_filt=detrend(Ch_data_filt);
Chan_Name={'Tx1Rx1 amp','TxRx1 ph','Tx2Rx1 amp','Tx2Rx1 ph','Tx3Rx1 amp','Tx3Rx1 ph','Tx4Rx1 amp','Tx4Rx1 ph',...
    'Tx1Rx2 amp','Tx1Rx2 ph','Tx2Rx2 amp','Tx2Rx2 ph','Tx3Rx2 amp','Tx3Rx2 ph','Tx4Rx2 amp','Tx4Rx2 ph',...
    'Tx1Rx3 amp','TxRx3 ph','Tx2Rx3 amp','Tx2Rx3 ph','Tx3Rx3 amp','Tx3Rx3 ph','Tx4Rx3 amp','Tx4Rx3 ph',...
    'Tx1Rx4 amp','TxRx4 ph','Tx2Rx4 amp','Tx2Rx4 ph','Tx3Rx4 amp','Tx3Rx4 ph','Tx4Rx4 amp','Tx4Rx4 ph','EOG1'}';
opts3.tWinBR = 12; % Window on which br is estimated

opts3.minInterceptDist = 0.2; % minimum time (s) between two intercepts default 0.05
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.4;
opts3.tWinVar=15;opts3.tWinpp=15;
opts3.tWinPk=10; % Window for peak detection moving average

% VarFeature=[VarBR_mean; VarPP_mean; VarIn_mean; VarEx_mean; meanBR;  meanIn; meanEx];
% CoV = std/mean  *100   Cov is calculated in a time window 
% % BR is average of all BRs over the time window, not pick only one BR
ind=[1,2,11,12,21,22,31,32];
% 8 slef channels , change direction - + 
% because BR is calculated based on  Max peak

% select the best channel by min cov
sign=[1,-1];
opts3.fig=0;
for i=1:length(ind)
    
[br_temp{1},pk_temp{1},~,VarFeature_temp{1},h_temp(1)] = brEstAvg(Ch_data_filt(:,ind(i))* sign(1),fsDS,opts3,Chan_Name(ind(i)),Motion);
[br_temp{2},pk_temp{2},~,VarFeature_temp{2},h_temp(2)] = brEstAvg(Ch_data_filt(:,ind(i))* sign(2),fsDS,opts3,Chan_Name(ind(i)),Motion);
    comp=[sum(VarFeature_temp{1}(1:2)),sum(VarFeature_temp{2}(1:2))];
    Opt=(find (comp==min(comp)));
    signOpt(i)= sign(Opt);
br{i}=br_temp(Opt);
pk{i}=pk_temp(Opt);
VarFeature{i}=VarFeature_temp(Opt);
h(i)=h_temp(Opt);
end


[br_eog,pk_eog,~,VarFeature_eog,h_eog] = brEstAvg(Ch_data_filt(:,33),fsDS,opts3,Chan_Name(33),Motion);



    % plot all including EOG
ind_all=[ind 33];
h_w=plotWaveform_s(Ch_data_filt(:,ind_all),Chan_Name(ind_all),fsDS);



%%  selection of channel 
for i=1:length(ind)
    temp=cell2mat(VarFeature{i});
    cov_comp(i)=sum(temp(1:2));

end
% sort the sequence of cov from the minium, the best channel is the first index
[cov_comp_seq,chOpt_seq] = sort(cov_comp); 

% use the first and second best ch 
Ch_opt1_num=chOpt_seq(1);
Ch_opt1_data=Ch_data_filt(:,ind(chOpt_seq(1)))*signOpt(chOpt_seq(1));
br_opt1_data=br(chOpt_seq(1));

Ch_opt2_num=chOpt_seq(2);
Ch_opt2_data=Ch_data_filt(:,ind(chOpt_seq(2)))*signOpt(chOpt_seq(2));
br_opt2_data=br(chOpt_seq(2));
%%

NCS_temp=Ch_opt1_data;
EMG_temp=Ch_data_filt(:,33);
[r1,p1] = corrcoef(normalize(NCS_temp),normalize(EMG_temp)); 
r_cor = r1(1,2);
[x,r]=xcorr(NCS_temp,EMG_temp,fsDS,'coeff');  % x+m, y 
max_coef=max(x);
max_r=r(find(x==max_coef))/fsDS;
h_cor=figure();
stem(x,r)
title(['Max Cor:',num2str(max_coef,'%.2f'),' TimeLag (s):',num2str(max_r,'%.3f')],'FontSize',11)

%%
status = mkdir([dataPath,'fig_REM']);

for i=1:2
figName = [dataPath,'fig_REM','\',fileName,'_Opt',num2str(i)];

print(h(chOpt_seq(i)),[figName,'.tiff'],'-dtiff','-r300');
savefig(h(chOpt_seq(i)),[figName,'.fig']);
end

figName = [dataPath,'fig_REM','\',fileName,'_waveform'];
print(h_w,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_w,[figName,'.fig']);

figName = [dataPath,'fig_REM','\',fileName,'_correlationPlot'];
print(h_cor,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_cor,[figName,'.fig']);
figName = [dataPath,'fig_REM','\',fileName,'_eog'];
print(h_eog,[figName,'.tiff'],'-dtiff','-r300');
savefig(h_eog,[figName,'.fig']);

%%

function [ampCh_filt_norm,ampCh_filt]=filtSignal(ampCh,toff,opt,fs,fsDS)
    
   
  
    ampCh=resample(ampCh,fsDS,fs);
    ampCh=ampCh((toff(1)*fsDS):toff(size(toff))*fsDS);
    ampCh_filt = filterLpHp(ampCh,fsDS,opt); % th amp
    
    ampCh_filt_norm = normalize(ampCh_filt);

end 

function [ampCh_filt_norm,ampCh_filt,ampCh_en,ampCh_raw]=filtSignal_emg(ampCh,toff,opt)
    
    fs=1e3;  % biopac sampling rate =1e3 
    fsDS=opt.fsDS; % down sampling rate  same as NCS  
    
  
    ampCh=resample(ampCh,fsDS,fs);
    ampCh_raw=ampCh((toff(1)*fsDS):toff(size(toff))*fsDS);
    
    [ampCh_en,~] = envelope(ampCh_raw,5,'peak');
    
    ampCh_en=detrend(ampCh_en);
    ampCh_filt = filterLpHp(ampCh_en,fsDS,opt); % th amp
    
    
    ampCh_filt_sm = sgolayfilt(ampCh_filt,opt.SMorder,opt.SMlen);
    %ampCh_filt_sm = smoothdata(ampCh_filt,'movmean',200);
    
    ampCh_filt_norm = normalize(ampCh_filt_sm);
    
end 

