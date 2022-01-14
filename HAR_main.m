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
arffwrite_HAR('HAR_52participants', AllFeatures, Name)
