function [CM, freq]=retrieve_RMap_GOI(cond, conn_metric, conn_nodes, dataset, frequencyrange, ROIs, norm, mode, GOI_selected, atlas, CMavg, CM)

% GOI = 'string' with gene of interest symbol name
% atlas = 'Desikan-Killiany' or
% condition = 'healthy' or 'PD'
% connectivity_metric = 'coh','icoh','grang'
% connectivitynodes = 186 or 1028
% dataset = 'OMEGA' or 'HCP'
% frequencyrange = 'all' (3:87 or 3:36), [f1 : f2], or 'theta', 'alpha', 'beta low','beta high','beta total','gamma narrow'
% ROIs = 'all', {'M1r','M1l','SMA','preSMA','IFGr','IFGl','mediotemporal left'; 'occipital pole left'} or [X Y Z]
% norm = '_norm' or '' %spatially normalized data
% mode = 'whole-brain connectivity', 'spectral plot', 'matrix'
addpath('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\3_TUTTI\rawdataCM')
addpath('C:\Users\Jonathan\Documents\MATLAB\add_on_Matlab\spm12')
addpath(genpath('C:\Users\Jonathan\Documents\MATLAB\add_on_Matlab\lead'))
addpath('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\9_Toolboxes\jvh_toolbox')
addpath(genpath('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\9_Toolboxes\jvh_toolbox'))
addpath('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\9_Toolboxes\jvh_toolbox')

filename = ['CMavg',num2str(conn_nodes),'_',dataset,'_',cond,'_',conn_metric,norm];
load(filename,'frequencies','coorddouble')
disp(filename)
outputname = [filename '_' frequencyrange];


%need to fix this later as bypass. Grang should have correct freq stored, not coh freq
if strcmp(conn_metric,'grang')
    frequencies=frequencies(1:36);
    connectivity_metric_label = 'net Granger causality';
elseif strcmp(conn_metric,'coh')
    connectivity_metric_label = 'coherence';
elseif strcmp(conn_metric,'icoh')
    connectivity_metric_label = 'imaginary part of coherence';
else
    error('unknown connectivity metric')
end

if strcmp(norm,'_norm')
    error('genetic analysis not possible for normative data')
    norm_label = 'spatially normalized';
else
    norm_label = 'absolute averages';
end

% try
%     CM = CMavg; % nodes x nodes x freq
% catch
%     CM = CMavg_norm; % nodes x freq
% end

freq = determine_frequencies(frequencyrange, frequencies);
ROI  = ROIs.determine_ROI();
%% preset CM

% if strcmp(norm,'_norm')
%     error('genetic analysis not possible for normative data')
%     CM=squeeze(CM(ROI.index,:)); %CM is nodes x freq
% else
%     CM=squeeze(mean(CM(ROI.index,:,:),1)); %CM is nodes x freq
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\6_Dev\');
cd('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\6_Dev\');

%% read in gene of interest
gene=GOI_selected;

%% read in the gene expression map for Desikan
gene_expression_map_file = 'abagen_Desikan_expression_map.csv';

% the genetic map is so large, that csvread does not work. Use datastore
ds = datastore(gene_expression_map_file, 'MissingValue',0); 

% retrieve expression for GOI
ds.SelectedVariableNames = gene;
gene_expression_map = readall(ds);
gene_expression_map = table2array(gene_expression_map);

%% read in how to project the CM to the Desikan
if conn_nodes==186
    load('channel_2_parcelnumber_Desikan_186chans.mat')
    channel_2_parcelnumber;
else
    error('not implemented for 1028 nodes')
end

%% load atlas
if strcmp(atlas,'Desikan-Killiany')
    nii_atlas = ea_load_nii('Desikan_cortical.nii');
end
% %% read in the cortical CM
% load('CMavg186_OMEGA_healthy_coh.mat') % from C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\3_TUTTI\rawdataCM

%% CM rearranged in connectivity per parcel
CM = CM;%CM is nodes x FOI %CM(:,freq.frequency_indices{:});
connectivity_in_parcels = [] ;
for el=1:length(gene_expression_map(:,1))
    connectivity_in_parcels(el)=mean(CM(channel_2_parcelnumber==el)); %average across the frequencies
end

%% correlation
% R= coefficient
% P= p-values
% RL= lower bound interval  for each coefficient
% RU= upper bound interval  for each coefficient

%% plot figure of x(coherence) y(gene expression) with each dot = ROI
gene_expression_map_no_nan=gene_expression_map;
connectivity_in_parcels_no_nan=connectivity_in_parcels;
gene_expression_map_no_nan(isnan(connectivity_in_parcels))=[];
connectivity_in_parcels_no_nan(isnan(connectivity_in_parcels))=[];

[R,P,RL,RU] = corrcoef(connectivity_in_parcels_no_nan,gene_expression_map_no_nan,'Rows','complete');

%corrplot([connectivity_in_parcels,gene_expression_map'])
figure
scatter(connectivity_in_parcels_no_nan,gene_expression_map_no_nan')
xlabel({['averaged connectivity per node ', connectivity_metric_label], ' ', outputname}, 'interprete','none')
ylabel(['spatial gene expresssion of ', gene])
lsline
str=sprintf(' r = %1.2f \n p = %1.3f',R(1,2),P(1,2));
T = text(max(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'right');

%% compute and save images with RMAP and RPMAP (masked for p<0.05)
assert(length(freq.frequency_indices)==1, 'specify a specific range of frequencies')
CMavg=mean(CMavg(:,:,freq.frequency_indices{:}),3); %nodes x nodes x freq

connectivity_62x186 = zeros(62,186);
for el=1:length(gene_expression_map(:,1))
    connectivity_62x186(el,:,:) = mean(CMavg(:,channel_2_parcelnumber==el),2);
end
for el=1:length(gene_expression_map(:,1))
    connectivity_62x62(:,el,:) = mean(connectivity_62x186(:,channel_2_parcelnumber==el),2);
end

r_all = zeros(62,1);
p_all = zeros(62,1);
for el=1:length(gene_expression_map(:,1))
    [R,P,RL,RU] = corrcoef(connectivity_62x62(:,el),gene_expression_map, 'Row','complete');
    r_all(el)=R(1,2);
    p_all(el)=P(1,2);
end
nii_atlas.img = changem(nii_atlas.img,1000,0);
nii_atlas.img = changem(nii_atlas.img,r_all,[1:62]);
nii_atlas.img = changem(nii_atlas.img,0,NaN);
nii_atlas.img = changem(nii_atlas.img,NaN,1000);
nii_atlas.fname = [outputname, '_RMap_', gene, '.nii' ]; %'CMavg186_OMEGA_healthy_coh_Rmap_PARK7.nii';
ea_write_nii(nii_atlas);

% mask by p-values

r_all(p_all>0.05)=0;
nii_atlas = ea_load_nii('Desikan_cortical.nii');
nii_atlas.img = changem(nii_atlas.img,1000,0);
nii_atlas.img = changem(nii_atlas.img,r_all,[1:62]);
nii_atlas.img = changem(nii_atlas.img,0,NaN);
nii_atlas.img = changem(nii_atlas.img,NaN,1000);
nii_atlas.fname = [outputname, '_RPMap_', gene, '.nii' ]; %'CMavg186_OMEGA_healthy_coh_RPmap_PARK7.nii';
ea_write_nii(nii_atlas);

% get gene expression

nii_atlas = ea_load_nii('Desikan_cortical.nii');
nii_atlas.img = changem(nii_atlas.img,1000,0);
nii_atlas.img = changem(nii_atlas.img,gene_expression_map,[1:62]);
nii_atlas.img = changem(nii_atlas.img,0,NaN);
nii_atlas.img = changem(nii_atlas.img,NaN,1000);
nii_atlas.fname = ['Gene_expression_', gene, '.nii' ]; %'Gene_expression_PARK7.nii';
ea_write_nii(nii_atlas);

vis3d(nii_atlas.img,'viridis') %just continue after x tick error

end
