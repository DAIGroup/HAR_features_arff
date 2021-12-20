%% _____________MAIN CODE_____________%%
clc
close all

users_list = dir('dataset');
file_list=[];
for nUser=3:size(users_list)
    files=dir(['dataset/', users_list(nUser).name]);
    for iFile=1:size(files)
        if files(iFile).isdir==0
            file_list=[file_list;files(iFile)];    
        end
    end
end
% Function - Features Computation 
[WFeat, Name, Weka, AllFeatures] = HAR_Features (file_list); 

%%FILE .ARFF for Weka
arffwrite_HAR('HAR_Weka_Pau+Angela_52', AllFeatures, Name)

%% ________________PLOTS______________%%

%if (0)
      
%     assetime = 1/fs : 1/fs : ((1/fs)* length(mag_acc_m));
%     figure;
%     plot(assetime, mag_acc_m);
%     ylabel('Acceleration (m/s^2)');
%     xlabel('Samples');
%     
%     figure;
%     plot(acc_x,'r')
%     plot(assetime, acc_x)
%     hold on
%     plot(acc_y,'b')
%     plot(assetime, acc_y)
%     hold on
%     plot(acc_z, 'k');
%     plot(assetime, acc_z);
%     xlabel('Samples');
%     ylabel ('Acceleration (m/s^2)');
%     legend ('Acc_x','Acc_y','Acc_z')
%     end