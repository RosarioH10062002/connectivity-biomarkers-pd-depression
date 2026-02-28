clearvars;clc
close all;
%% CHECK THE PARTICIPANTS' EXISTANCE 
function check_subjects(subject_list, root_path, filename_expected)
    for k = 1:length(subject_list)
    
        subject = subject_list(k);
    
        fullpath = fullfile(root_path, subject, subject, ...
                            "data", filename_expected);
    
        if isfile(fullpath)
            fprintf("%s OK\n", subject);
        else
            fprintf("%s MISSING FILE\n", subject);
        end
    
    end
end 
%% ROI + GRANGER CAUSALITY 

function G_all = granger_causality(subject_list, root_path,filename_expected)
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
%% CHECK AND LOAD PARTICIPANTS CTRL 
clc; 
root_path = "CONN_PREPROCESSED_CTRL/16_participants";
status = "ctrl";
filename_expected = "ROI_Subject001_Session001.mat";
subject_list = [
    "101195", "102447", "103161", "103183", "103467", "103542", ...
    "108909", "113050", "115698", ...
    "130028", "138022", "142086", "149120", ...
    "149716", "156484", "160890"
];
check_subjects(subject_list,root_path, filename_expected); 
G_all = granger_causality(subject_list, root_path,filename_expected);
G_all_ctrl = G_all;

%% CHECK AND LOAD PARTICIPANTS PD 
clc; 
root_path = "CONN_PREPROCESSED_PD/16_participants";
subject_list = [
    "100005", "100006", "100007", "100268", "101018", ...
    "101146", "101174", "101476", "101479", "101735", "101742", ...
    "101751", "101841", "102012", "102053", "102078"
];
status = "PD"; 
filename_expected = "ROI_Subject001_Session001.mat";
check_subjects(subject_list,root_path, filename_expected);
G_all = granger_causality(subject_list, root_path,filename_expected);
G_all_pd = G_all;

%% CHECK AND LOAD PARTICIPANTS PD DEP 
clc;
root_path = "CONN_PREPROCESSED_PD_DEP/16_participants";
subject_list = [
    "100018", "100267", "100842", "101050", "102978", "133507", ...
    "139859", "142879", "144131", "149116", "172260", "174141", ...
    "182427", "184432", "211482", "219411"
];
status = "PD_DEP"; 
filename_expected = "ROI_Subject001_Session001.mat";
check_subjects(subject_list,root_path, filename_expected);
G_all = granger_causality(subject_list, root_path,filename_expected);
G_all_pdd = G_all;

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
    [AIC,BIC] = tsdata_to_infocrit(input, maxorder, 'LWR');
    
    %figure; plot(1:maxorder, AIC, '-o'); hold on;
    %plot(1:maxorder, BIC, '-o'); legend('AIC','BIC');
    %xlabel('Model order p'); ylabel('Information criterion');
   
    [~, p] = min(BIC);   % (more conservative for fMRI)
    p = 2; 
    disp(["Chosen order (BIC):", p]);

    [A,SIG] = tsdata_to_var(input, p, 'LWR');
    
    F = var_to_autocov(A, SIG);
    G = autocov_to_pwcgc(F);   % pairwise-conditional GC
  
end 

%% MULTIVARIATE GRANGER CAUSALITY FIGURES 
function mvgc_graph(G, roi_names_labels, subject, status, global_G, p)
    subject = string(subject); 
    status = string(status); 
    %tic
    f = figure('Visible','off');   
    
    imagesc(G)
    axis square
    clim([0 global_G])
    colorbar
    
    xticks(1:length(roi_names_labels))
    yticks(1:length(roi_names_labels))
    xticklabels(roi_names_labels)
    yticklabels(roi_names_labels)
    xtickangle(90)
    
    %title('MVGC (p=1)')
    title(sprintf('MVGC (p=%d) - Subject %s - %s', p, subject, status))
    filename = sprintf('MVGC_Subject_%s_p%d_%s.png', subject, p, status);
    exportgraphics(f, filename, 'Resolution', 300)
    
    close(f)  
    %toc 
end 
%% STATISTICAL ANALYSIS 
G_ctrl = cat(3, G_all_ctrl{:});
G_pd   = cat(3, G_all_pd{:});
G_pdd  = cat(3, G_all_pdd{:});

[nROI, ~, ~] = size(G_ctrl);
p_kw = zeros(nROI, nROI);   % Kruskal–Wallis p-values

for i = 1:nROI
    for j = 1:nROI
        
        if i == j
            p_kw(i,j) = NaN; 
            continue
        end
       
        x1 = squeeze(G_ctrl(i,j,:));
        x2 = squeeze(G_pd(i,j,:));
        x3 = squeeze(G_pdd(i,j,:));
        
        values = [x1; x2; x3];
        groups = [ones(size(x1));
                  2*ones(size(x2));
                  3*ones(size(x3))];
        
        p_kw(i,j) = kruskalwallis(values, groups, 'off');
        
    end
end

disp(p_kw)
significant_edges_whithout_correction = p_kw < 0.05;
figure;
imagesc(significant_edges_whithout_correction)
axis square
colormap(gray)
colorbar

xticks(1:length(roi_names_labels))
yticks(1:length(roi_names_labels))
xticklabels(roi_names_labels)
yticklabels(roi_names_labels)
xtickangle(90)

title("Significant edges (p < 0.05, uncorrected)")
   
%% 
pvec = p_kw(:);
valid = ~isnan(pvec);

pvec_fdr = nan(size(pvec));
pvec_fdr(valid) = mafdr(pvec(valid), 'BHFDR', true);

p_kw_fdr = reshape(pvec_fdr, nROI, nROI);
significant_edges = p_kw_fdr < 0.05;
figure;
imagesc(significant_edges)
axis square
colormap(gray)
colorbar

xticks(1:length(roi_names_labels))
yticks(1:length(roi_names_labels))
xticklabels(roi_names_labels)
yticklabels(roi_names_labels)
xtickangle(90)

title("Significant edges (p < 0.05, corrected)")