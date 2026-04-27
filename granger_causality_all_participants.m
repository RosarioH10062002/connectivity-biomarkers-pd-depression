this is my actual code: clc;
root_path_ctrl = "G:\Mon Drive\M2\Impact Scholar Programme\Preprocessed_data\CONN_PREPROCESS_CTRL";
root_path_pd = "G:\Mon Drive\M2\Impact Scholar Programme\Preprocessed_data\CONN_PREPROCESS_PD";
root_path_pdd = "G:\Mon Drive\M2\Impact Scholar Programme\Preprocessed_data\CONN_PREPROCESS_PD_DEP";

list_part_pdd = clean_list_participants(root_path_pdd);
list_part_pd = clean_list_participants(root_path_pd);
list_part_ctrl = clean_list_participants(root_path_ctrl);

filename_expected = "ROI_Subject001_Session001.mat";

%% OBTAIN DATA 
%CTRL
G_ctrl = granger_causality_function(list_part_ctrl, root_path_ctrl,filename_expected);
G_pd = granger_causality_function(list_part_pd, root_path_pd,filename_expected);
G_pdd = granger_causality_function(list_part_pdd, root_path_pdd,filename_expected);

%% CLEAN PARTICIPANTS

function subject_list = clean_list_participants(root_path)
    files = dir(root_path);
    names = string({files.name});
    isSubject = ~isnan(str2double(names)); %str2double = string to number 
    subject_list = names(isSubject);
end

%% COMPUTE ROI 
function [final_data, roi_names_labels] = compute_roi(data_load)
    roi_names= string(data_load.names);
    idx_DMN = startsWith(roi_names,"networks.DefaultMode.");
    idx_SAL = startsWith(roi_names, "networks.Salience.");
    idx_FPN = startsWith(roi_names, "networks.FrontoParietal.");
    mask_regions = (idx_DMN|idx_SAL|idx_FPN); 
    roi = data_load.data(:,mask_regions); 
    roi_names_labels = roi_names(:,mask_regions);
    roi_names_labels = replace(roi_names_labels, "networks.DefaultMode.", "DMN-");
    roi_names_labels = replace(roi_names_labels, "networks.Salience.", "SAL-"); 
    roi_names_labels = replace(roi_names_labels, "networks.FrontoParietal.", "FPN-");
    final_data = cell2mat(roi);
end 
%% MULTIVARIATE GRANGER CAUSALITY
close all;
addpath(genpath('path_to_mvgc'))
startup

function [G,p] = mvgc(final_data)
    input = final_data'; 
    maxorder = 5; 
    [~,BIC] = tsdata_to_infocrit(input, maxorder, 'LWR');
    
    %figure; plot(1:maxorder, AIC, '-o'); hold on;
    %plot(1:maxorder, BIC, '-o'); legend('AIC','BIC');
    %xlabel('Model order p'); ylabel('Information criterion');
   
    [~, p] = min(BIC);   % (more conservative for fMRI)
    disp(["Chosen order (BIC):", p]);

    [A,SIG] = tsdata_to_var(input, p, 'LWR');
    
    F = var_to_autocov(A, SIG);
    G = autocov_to_pwcgc(F);   % pairwise-conditional GC
  
end 
%% ROI + GRANGER CAUSALITY 

function G_all = granger_causality_function(subject_list, root_path,filename_expected)
    nSub = length(subject_list);
    
    G_all = cell(nSub,1); 
    P_all = zeros(nSub,1);
    subjects = strings(nSub,1);
    roi_labels_all = cell(nSub,1);
    max_G_all = zeros(nSub,1); 
    
    for k = 1:nSub
    
        subject = subject_list(k);
    
        fullpath = fullfile(root_path, subject, subject, ...
                            "data", filename_expected);
     
        fprintf("Processing subject %s\n", subject);
    
        data_load = load(fullpath);
    
        [final_data, roi_names_labels] = compute_roi(data_load);
        [G,p] = mvgc(final_data); 
        %p = 2; 
    
        G_all{k} = G; 
        P_all(k) = p;
        subjects(k) = subject;
        max_G_all(k) = max(G(:));
        roi_labels_all{k} = roi_names_labels;
    
    end
    %{
    global_G = max(max_G_all);
    
    for k = 1:length(G_all)
    
        mvgc_graph(G_all{k}, roi_labels_all{k}, subjects(k), status, global_G, P_all(k));
    end
    %}
end 

%% SAVE MY DATA 
save('GC_results.mat', 'G_ctrl', 'G_pd', 'G_pdd')
