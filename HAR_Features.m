function [WFeat, Name, Weka, AllFeatures] = HAR_Features (file_list)
%% F E A T U R E S %%
Features = struct;
AllFeatures=[];
for iFile=1:length(file_list)
    filename = (file_list(iFile).name); % Import file (almost 5 minutes of acquisition)
    disp(filename);
    data = csvread([file_list(iFile).folder, '/', file_list(iFile).name]); 
    data = data./64; 
    duration = length(data);

    %---Activity
    format='%u';
    STR=filename; %taken from the name of the file
    if contains(STR,'drink_water')  
        s(1)=1;
    elseif contains(STR,'eat_meal') 
        s(1)=2;
    elseif contains(STR,'open_a_bottle') 
        s(1)=3;
    elseif contains(STR,'open_a_box')
        s(1)=4;
    elseif contains(STR,'brush_teeth') 
        s(1)=5;
    elseif contains(STR,'brush_hair')
        s(1)=6;
    elseif contains(STR,'take_off_a_jacket')
        s(1)=7;
    elseif contains(STR,'put_on_a_jacket')
        s(1)=8;
    elseif contains(STR,'put_on_a_shoe')
        s(1)=9;
    elseif contains(STR,'take_off_a_shoe')
        s(1)=10;
    elseif contains(STR,'put_on_glasses')
        s(1)=11;
    elseif contains(STR,'take_off_glasses')
        s(1)=12;
    elseif contains(STR,'sit_down')
        s(1)=13;
    elseif contains(STR,'stand_up')
        s(1)=14;
    elseif contains(STR,'writing')
        s(1)=15;
    elseif contains(STR,'phone_call')
        s(1)=16;
    elseif contains(STR,'type_on_a_keyboard')
        s(1)=17;
    elseif contains(STR,'salute')
        s(1)=18;
    elseif contains(STR,'sneeze_cough')
        s(1)=19;
    elseif contains(STR,'blow_nose')
        s(1)=20;
    elseif contains(STR,'washing_hands')
        s(1)=21;
    elseif contains(STR,'dusting')
        s(1)=22;
    elseif contains(STR,'ironing')
        s(1)=23;
    elseif contains(STR,'washing_dishes')
        s(1)=24;
    else
        display([filename, ' ', STR]);
        s(1)=99;  %Error handle if doesnt contain any of the previous
    end

    chunks=split(STR,'.');
    chunks=split(chunks{1,1},'_');
    sizeChunks=size(chunks);
    s(2)=str2num(chunks{sizeChunks(1)-1,1});
    s(3)=str2num(chunks{sizeChunks(1),1})+1;
    data = data - mean(data);  
    L = length(data(:,1)); 

    %% Filtering 
    %----BandPass Filter Butterworth - 4 order
    fs= 32; %sampling frequency
    f_cut= 15; %cut-off freq.
    [b,a] = butter(4,15/(fs/2),'low');
    acc = filter(b,a,data);

    %----Median filter - 3 order 
    acc_M = medfilt1(acc,3);
    acc_x = acc_M(:,1);     
    acc_y = acc_M(:,2);
    acc_z = acc_M(:,3);
    acc_m = [acc_x acc_y acc_z];

    %% ----Features in Time Domain
    %w = 3; %window in seconds
    %delta = round( w*fs ); %window in samples
    time_w = min(5*fs,duration); %windows in seconds
    overlap = min(1*fs,duration); %overlap in seconds
    mag_time = sqrt(acc_x.^2 + acc_y.^2 + acc_z.^2); %Signal Magnitud Vector (MAGNITUDE)

    clear Mean;
    clear Median;
    clear Standev;
    clear Maximum;
    clear Minimum;
    clear Range;
    clear Corraxis;
    clear SMA;
    clear CV;
    clear MedAbsDev;
    clear Skewness ;
    clear Kurtosis ;
    clear Autocorr;
    clear Pt_20;
    clear Pt_50;
    clear Pt_80;
    clear IQ;
    clear N_Peaks;
    clear Peak2Peak_Amplitude;
    clear Energy;
    clear Pt_90;
    clear RMS;
    % Feature in Frequency
    clear Entropy ;
    clear SEnergy;
    clear SCentroid;
    clear MeanF;
    clear StdevF;
    clear Pt_25F;
    clear Pt_50F;
    clear Pt_75F;
        
    i=0;
    for p1=1:overlap:duration
        p2= p1 + time_w-1;
        if p2>duration
            break
        end
        i=i+1;
        Mean(i,1)= mean(acc_x(p1:p2)); %x
        Mean(i,2)= mean(acc_y(p1:p2)); %y
        Mean(i,3)= mean(acc_z(p1:p2)); %z
        Mean(i,4)= mean(mag_time(p1:p2)); %MAG

        Median(i,1)= median(acc_x(p1:p2));
        Median(i,2)= median(acc_y(p1:p2));
        Median(i,3)= median(acc_z(p1:p2));
        Median(i,4)= median(mag_time(p1:p2));

        Standev(i,1)= std(acc_x(p1:p2));
        Standev(i,2)= std(acc_y(p1:p2));
        Standev(i,3)= std(acc_z(p1:p2));
        Standev(i,4)= std(mag_time(p1:p2));

        Maximum(i,1)= max(acc_x(p1:p2));
        Maximum(i,2)= max(acc_y(p1:p2));
        Maximum(i,3)= max(acc_z(p1:p2));
        Maximum(i,4)= max(mag_time(p1:p2));

        Minimum(i,1)= min(acc_x(p1:p2));
        Minimum(i,2)= min(acc_y(p1:p2));
        Minimum(i,3)= min(acc_z(p1:p2));
        Minimum(i,4)= min(mag_time(p1:p2));

        Range(i,1)= range(acc_x(p1:p2));
        Range(i,2)= range(acc_y(p1:p2));
        Range(i,3)= range(acc_z(p1:p2));
        Range(i,4)= range(mag_time(p1:p2));

        Corraxis(i,1)= corr(acc_x(p1:p2),acc_y(p1:p2));
        Corraxis(i,2)= corr(acc_y(p1:p2),acc_z(p1:p2));
        Corraxis(i,3)= corr(acc_z(p1:p2),acc_x(p1:p2));

        %Signal Magnitude Area
        SMA(i,1)= sum(abs(acc_x(p1:p2)) + abs(acc_y(p1:p2)) + abs(acc_z(p1:p2)));

        %Coefficient of Variation
        CV(i,1)= std(acc_x(p1:p2))/Mean(i,1);
        CV(i,2)= std(acc_y(p1:p2))/Mean(i,2);
        CV(i,3)= std(acc_z(p1:p2))/Mean(i,3);
        CV(i,4)= std(mag_time(p1:p2))/Mean(i,4);

        %Median Absolute Deviation
        MedAbsDev(i,1)= sum(abs(acc_x(p1:p2)-Median(i,1)))/ L;
        MedAbsDev(i,2)= sum(abs(acc_y(p1:p2)-Median(i,2)))/ L;
        MedAbsDev(i,3)= sum(abs(acc_z(p1:p2)-Median(i,3)))/ L;
        MedAbsDev(i,4)= sum(abs(mag_time(p1:p2)-Median(i,4)))/length(mag_time);

        %Skewness
        Skewness(i,1)= skewness(mag_time(p1:p2));
        % Kurtosis
        Kurtosis(i,1)= kurtosis(mag_time(p1:p2));

        %Autocorrelation
        for k= 1:length(mag_time(p1:p2))-1
            Autocorr(i,1)= (sum((mag_time(k)-Mean(i,4))*(mag_time(k+1)-Mean(i,4))))/sum((mag_time(k)-Mean(i,4))^2);
        end
    
        %Percentiles (20,50,80) and Interquartile Range
        Pt_20(i,1)= prctile(mag_time(p1:p2),20);
        Pt_50(i,1)= prctile(mag_time(p1:p2),50);
        Pt_80(i,1)= prctile(mag_time(p1:p2),80);
        IQ(i,1)= prctile(mag_time(p1:p2),75)- prctile(mag_time(p1:p2),25);

        %Location Peak
        [peaks,locs,w]= findpeaks(mag_time(p1:p2)); %returns peaks, index e width
        N_Peaks(i,1)= length(peaks);
        Peak2Peak_Amplitude(i,1)= max(w)-min(w);

        %Energy 
        Energy(i,1) = sum((mag_time(p1:p2)).^2)/(fs/L);
        %Mozos' features
        Pt_90(i,1)= prctile(mag_time(p1:p2),90); %percentile 90
        % Pt_90(i,2)= prctile(acc_x(p1+1:p2),90); 
        % Pt_90(i,3)= prctile(acc_y(p1+1:p2),90);
        % Pt_90(i,4)= prctile(acc_z(p1+1:p2),90); 
        RMS(i,1)= rms(mag_time(p1:p2)); %root mean square
        %deltax(i,1)= (1/(length(L1)-1))*sum(abs(mag_time(p1+1:p2)));%means of the abs values of first differences of raw signal
        %the means of the absolute values of the second differences of the raw signals
    
        %% ----Frequency Domain
        xFFT = fft(acc_x(p1:p2));
        yFFT = fft(acc_y(p1:p2));
        zFFT = fft(acc_z(p1:p2));
        magFFT= fft(mag_time(p1:p2));
        %magFFT_half = magFFT(1:round(length(magFFT))/2);
        magAbs= abs(magFFT);
        %magPhase = (angle(magFFT)).*57.295779513;
        
        %% ----Features in Frequency Domain
        %Spectral Entropy
        Entropy (i,1) = wentropy(magAbs,'shannon');
        %Spectral Energy
        SEnergy (i,1) = sum((magAbs).^2)/time_w;
        %Spectral Centroid
        SCentroid (i,1) = mean(magAbs)/time_w; 

        %Mozos' features in freq
        pxx = periodogram(magAbs);
        px = periodogram(xFFT);
        py = periodogram(yFFT);
        pz = periodogram(zFFT);
        MeanF(i,1)= mean(pxx);
        MeanF(i,2)= mean(px);
        MeanF(i,3)= mean(py);
        MeanF(i,4)= mean(pz);
        StdevF(i,1)= std(pxx);
        StdevF(i,2)= std(px);
        StdevF(i,3)= std(py);
        StdevF(i,4)= std(pz);
        Pt_25F(i,1)= real(prctile(magAbs,25));
        % Pt_25F(n,2)= real(prctile(xFFT(t1+1:t2),25)); %percentile 25
        % Pt_25F(n,3)= real(prctile(yFFT(t1+1:t2),25)); 
        % Pt_25F(n,4)= real(prctile(zFFT(t1+1:t2),25)); 
        Pt_50F(i,1)= real(prctile(magAbs,50)); %percentile 50
        % Pt_50F(n,2)= real(prctile(xFFT(t1+1:t2),50));
        % Pt_50F(n,3)= real(prctile(yFFT(t1+1:t2),50));
        % Pt_50F(n,4)= real(prctile(zFFT(t1+1:t2),50));
        Pt_75F(i,1)= real(prctile(magAbs,75));  %percentile 75
        % Pt_75F(n,2)= real(prctile(xFFT(t1+1:t2),75));
        % Pt_75F(n,3)= real(prctile(yFFT(t1+1:t2),75));
        % Pt_75F(n,4)= real(prctile(zFFT(t1+1:t2),75)); 
    end
    
    %% STRUCTURE FOR FEATURES

    WFeat.Mean = Mean;
    WFeat.Median = Median;
    WFeat.Stdev = Standev;
    WFeat.Max = Maximum;
    WFeat.Min = Minimum;
    WFeat.Range = Range;
    WFeat.Corraxis = Corraxis;
    WFeat.Sma = SMA;
    WFeat.Cv = CV;
    WFeat.MedAbsDev = MedAbsDev;
    WFeat.Skewness = Skewness ;
    WFeat.Kurtosis = Kurtosis ;
    WFeat.Autocorr = Autocorr;
    WFeat.Pt20 = Pt_20;
    WFeat.Pt50 = Pt_50;
    WFeat.Pt80 = Pt_80;
    WFeat.Iq = IQ;
    WFeat.NPeaks = N_Peaks;
    WFeat.Peak2PeakAmp = Peak2Peak_Amplitude;
    WFeat.Energy = Energy;
    WFeat.Pt_90=Pt_90;
    WFeat.RMS=RMS;
    % Feature in Frequency
    WFeat.Entropy = Entropy ;
    WFeat.SEnergy = SEnergy;
    WFeat.Centroid = SCentroid;
    WFeat.MeanF = MeanF;
    WFeat.StdevF = StdevF;
    WFeat.Pt_25F= Pt_25F;
    WFeat.Pt_50F=Pt_50F;
    WFeat.Pt_75F=Pt_75F;
    Features.data = WFeat;

    %% MATRIX FOR FEATURES in WINDOWS
    MatrixFeat_W = [ Mean, Median, Standev, Maximum,...
               Minimum, Range, Corraxis, SMA...
               CV, MedAbsDev, Skewness, Kurtosis,...
               Autocorr, Pt_20, Pt_50, Pt_80, IQ,...
               N_Peaks, Peak2Peak_Amplitude, Energy, Pt_90, RMS, Entropy ,...
               SEnergy, SCentroid, MeanF, StdevF, Pt_25F, Pt_50F,Pt_75F];

    Name = {'MeanX', 'MeanY' , 'MeanZ', 'MeanMag', 'MedianX', 'MedianY', 'MedianZ',...
        'MedianMag', 'StdevX','StdevY', 'StdevZ', 'StdevMag', 'MaximumX', 'MaximunY', 'MaximumZ', 'MaximumMag',...
        'MinimumX', 'MinimumY', 'MinimumZ', 'MinimumMag', 'RangeX', 'RangeY', 'RangeZ', 'RangeMag',...
        'CorraxisXY','CorraxisYZ', 'CorraxisZX','SMA', 'CvX', 'CvY', 'CvZ', 'CvMag',...
        'MedAbsDevX', 'MedAbsDevY', 'MedAbsDevZ', 'MedAbsDevMag','SkewnessMag', 'KurtosisMag', ...
        'AutocorrMag','Pt20',  'Pt50', 'Pt80', 'Iq', 'NPeaks','Peak2PeakAmp','Energy'...
        'Pt_90', 'RMS', 'Entropy', ...
        'SEnergy', 'SCentroid', 'MeanF', 'MeanFX', 'MeanFY' , 'MeanFZ',...
        'StdevF', 'StdevFX','StdevFY', 'StdevFZ','Pt_25F', 'Pt_50F','Pt_75F' }';

    Wfeatures = MatrixFeat_W;
    % Matrix for .arff file
    feat=repmat([s(1),s(2),s(3)],size(Wfeatures,1),1); %initial features:ID,activity,scenario, check
    Weka=[Wfeatures,feat];
    AllFeatures = [AllFeatures; Weka];

end
%assetime = 1/fs : 1/fs : ((1/fs)* length(data));
%     figure;
%     plot(assetime, mag_acc_m);
%     ylabel('Acceleration (m/s^2)');
%     xlabel('Samples');
    
%     figure;
%     %plot(acc_x,'r')
%     plot(assetime, data(1:9000,1))
%     hold on
%     %plot(acc_y,'b')
%     plot(assetime, data(1:9000,2))
%     hold on
%     %plot(acc_z, 'k');
%     plot(assetime, data(1:9000,3));
%     xlabel('Time [s]');
%     ylabel ('Acceleration [m/s^2]');
%     legend ('Acc_x','Acc_y','Acc_z')
%     %xlim([0 290])
