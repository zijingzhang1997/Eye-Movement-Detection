% feature=[mean(br), std(br), mean(pp), std(pp), mean(in), std(in), mean(ex), std(ex),... 
%          skew_mean, kurt_mean, entro,...
%            per_power, cycle, covBR, covPP, void_t, ...
%            max_in, max_ex, max_br, min_br, min_pp, ...
%            Cor_br, SD_br, Cor_pp, SD_pp, Cor_in, SD_in, Cor_ex, SD_ex];



function [featAll,pkidx]=FeatInEpoch(EpochData,fs,opt)
featAll=zeros(opt.featNum,size(EpochData,2)); 

for i=1:size(EpochData,2)
    EpochData_row=EpochData(:,i);
    
    [featAll(:,i),pkidx{i}]=featExtract(EpochData_row,fs,opt);
    
end  

end

    

function [feature,pkidx]=featExtract(Data,fs,opt)    
    

%feature=zeros(opt.featNum,1);   
t = (0:(length(Data)-1))/fs;
Data=detrend(Data);
Data = normalize(Data);
% Minima and maxima detection 



pk = findMaxMin(Data,fs,opt);
if isempty(pk)
    feature=zeros(opt.featNum,1); 
    feature(:)=1000;
    pkidx=[];
     return
end
if pk(1).ind(1) == 1
    % If first peak is inhalation peak, skip it
    pk(1).ind = pk(1).ind(2:end); 
    pk(1).idx = pk(1).idx(2:end);
end

if pk(1).ind(end) == 0
    pk(1).ind = pk(1).ind(1:end-1);
    pk(1).idx = pk(1).idx(1:end-1);
end
pkMax = pk(1).idx(pk(1).ind == 1);
pkMin = pk(1).idx(pk(1).ind == 0);

%% calibrate peaks 
if opt.calib==1
  
   
    del = zeros(length(pkMax),1);

   
    for i = 1:length(pkMax)
            del(i) = abs(Data(pkMax(i))-Data(pkMin(i))); % Making it positive always
    end

    tBRmax1 = t(pkMax);
    idxCalib = ((tBRmax1 >= opt.calibT(1))&(tBRmax1 <= opt.calibT(2)));
    delNcsCalib = mean(del(idxCalib)); 


    idxPk = (del >= opt.calibMinPkRatio*delNcsCalib);
    nold=length(pkMax);
    pkMax = pkMax(idxPk); % Only keep indices that satisfy, otherwise ignore
    pkMin = pkMin(idxPk);
    nn=length(pkMax);
    CaliDel=nold-nn;
    
end

    %% Find features in each cycle (min-max-min)
cycle=length(pkMax)-1;
br = zeros(cycle,1);ibi = zeros(cycle,1);pp = zeros(cycle,1);in = zeros(cycle,1);ex = zeros(cycle,1);
ibi = zeros(cycle,1);  % inter breath time interval 
IEpp= zeros(cycle,1);  % inhale volume / exhale volume 
IER= zeros(cycle,1);  % inhale time / exhale time ratioe 
skew=zeros(cycle,1);kurt=zeros(cycle,1);
for i=1:cycle
        br(i)=60/(t(pkMin(i+1))-t(pkMin(i)));  %BR unit BPM 
        ibi(i)=t(pkMin(i+1))-t(pkMin(i));  %Inter-breath interval
        pp(i)=(Data(pkMax(i))-Data(pkMin(i))-Data(pkMin(i+1)))/2;  %% peak to peak averaged on min-max & max-min 
        IEpp(i) =(Data(pkMax(i))-Data(pkMin(i)))/(Data(pkMax(i))-Data(pkMin(i+1)));
    
        in(i)=t(pkMax(i))-t(pkMin(i));
        ex(i)=t(pkMin(i+1))-t(pkMax(i));
        IER(i)=in(i)/ex(i);
        
        
        DataSeg=Data(pkMin(i):pkMin(i+1));
        kurt(i)=kurtosis(DataSeg);  %measure of the "tailedness"
        skew(i)=skewness(DataSeg);  %measure of the asymmetry
        en(i)=entropy(DataSeg);
    end
 
entro=entropy(Data); 
skew_mean=mean(skew);kurt_mean=mean(kurt);

covBR=std(br)/mean(br);
covIBI=std(ibi)/mean(ibi);
covPP=std(pp)/mean(pp);
covIN=std(in)/mean(in);
covEX=std(ex)/mean(ex);

% auto correlation and successive difference in neighbor cycles 
Cor_br = xcorr(br,1,'coeff');Cor_br =Cor_br (1);
SD_br= mean(abs(diff(br))./br(1:end-1));
Cor_ibi = xcorr(ibi,1,'coeff');Cor_ibi =Cor_ibi (1);
SD_ibi= mean(abs(diff(ibi))./ibi(1:end-1));
Cor_pp = xcorr(pp,1,'coeff');Cor_pp =Cor_pp (1);
SD_pp= mean(abs(diff(pp))./pp(1:end-1));
Cor_in = xcorr(in,1,'coeff');Cor_in =Cor_in (1);
SD_in= mean(abs(diff(in))./in(1:end-1));
Cor_ex = xcorr(ex,1,'coeff');Cor_ex =Cor_ex (1);
SD_ex= mean(abs(diff(ex))./ex(1:end-1));
Cor_ex = xcorr(ex,1,'coeff');Cor_ex =Cor_ex (1);
SD_ex= mean(abs(diff(ex))./ex(1:end-1));
Cor_IEpp = xcorr(IEpp,1,'coeff');Cor_IEpp =Cor_IEpp (1);
SD_IEpp= mean(abs(diff(IEpp))./IEpp(1:end-1));
Cor_IER = xcorr(IER,1,'coeff');Cor_IER =Cor_IER (1);
SD_IER= mean(abs(diff(IER))./IER(1:end-1));
     
%% frequency features 
% SNR br _hr_
if ~isempty(opt.patient_rr)
patient_rr=opt.patient_rr;  % use mean BR from feature extraction 
patient_hr=opt.patient_hr; % use reference HR
else 
patient_rr=18;  % use mean BR from feature extraction 
patient_hr=80; % use reference HR
end
L = length(Data);
Y = fft((Data));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
noise= median(P1);

signal_rr_max = max(P1(f>(patient_rr/60-0.2) & f<(patient_rr/60+0.2)));
signal_hr_max = max(P1(f>(patient_hr/60-0.2) & f<(patient_hr/60+0.2)));
signal_rr=f(P1==signal_rr_max);  %   frequency with maximum energy for RR ; need to change to BPM
signal_hr=f(P1==signal_hr_max);  % frequency with maximum energy for HR ; need to change to BPM
snr_hr = 10*log10(signal_hr_max/noise);
snr_br = 10*log10(signal_rr_max/noise);
%% spectrogtam 
[spec_plot,f_spec,t_spec] = spectrogram(detrend(Data/noise), kaiser(fs*20,5), fs*5, 2^nextpow2(fs*20), fs, 'yaxis');

freq=opt.freq_range;   % frequency for different range of spec. 
%Time averaged power density (dB/Hz):
power_spec(1) = mean(10*log10(abs(spec_plot(f_spec>freq(1,1) & f_spec<freq(1,2)))),'all');
power_spec(2) = mean(10*log10(abs(spec_plot(f_spec>freq(2,1) & f_spec<freq(2,2)))),'all');
power_spec(3) = mean(10*log10(abs(spec_plot(f_spec>freq(3,1) & f_spec<freq(3,2)))),'all');

fr=opt.freq_range_rr;  % calculate power density in frequency range [-fr,+fr]close to reported HR, RR from PSG 
power_spec(4) = mean(10*log10(abs(spec_plot(f_spec>(signal_rr-fr) & f_spec<(signal_rr+fr)))),'all');
power_spec(5) = mean(10*log10(abs(spec_plot(f_spec>(signal_hr-fr) & f_spec<(signal_hr+fr)))),'all');

% power ratio in frequency range /all power 
[pxx,f] = periodogram(Data,hamming(length(Data)),length(Data),fs,'power');


pTot = bandpower(pxx,f,'psd');
power_ratio(1) = 100*(bandpower(pxx,f,[freq(1,1) freq(1,2)],'psd')/pTot);
power_ratio(2) = 100*(bandpower(pxx,f,[freq(2,1) freq(2,2)],'psd')/pTot);
power_ratio(3) = 100*(bandpower(pxx,f,[freq(3,1) freq(3,2)],'psd')/pTot);
power_ratio(4) = 100*(bandpower(pxx,f,[max(0,signal_rr-fr) signal_rr+fr],'psd')/pTot);
power_ratio(5) = 100*(bandpower(pxx,f,[max(0,signal_hr-fr) signal_hr+fr],'psd')/pTot);

%% 

feature_t=[mean(br),std(br),mean(ibi),std(ibi),mean(pp),std(pp),mean(in),std(in),mean(ex),std(ex),...
    mean(IEpp),std(IEpp),mean(IER),std(IER),...
     ];
        
       
feature=[feature_t  feature_f];     
feature(isnan(feature))=100;  % NA  change to 100, means error is very big 



pkidx.max=pkMax;    
pkidx.min=pkMin; 



% if using feature detection can't work, error!  Then, directly give value 



    
end



