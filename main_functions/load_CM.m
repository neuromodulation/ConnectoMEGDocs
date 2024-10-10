function [ba, freq, ROIs, CMfreq, CMroi, CMavg, filename] = load_CM(condition, connectivity_metric, conn_nodes, significant, dataset, frequencyrange, ROIs, norm, previous_filename, previous_CMavg, previous_CMfreq, previous_CMroi, smoothening)
    arguments (Input)
        condition           {mustBeMember(condition,['healthy','PD'])} ='healthy'
        connectivity_metric {mustBeMember(connectivity_metric,['coh','icoh','grang'])} = 'coh'
        conn_nodes          {mustBeMember(conn_nodes,[186,1018,1028])} = 186
        significant         {mustBeA(significant, "logical"), mustBeScalarOrEmpty} = false
        dataset             {mustBeMember(dataset,['OMEGA','HCP'])} = 'OMEGA'
        frequencyrange      {mustBeValidFrequencyrange} = 'all'
        ROIs                {mustBeA(ROIs,"ROI")} = struct();
        norm                {mustBeMember(norm,['_norm',''])} = ''
        previous_filename   {mustBeText} = ''
        previous_CMavg      double = []
        previous_CMfreq     double = []
        previous_CMroi      double = []
        smoothening         {mustBeInteger} = 0
    end
    arguments (Output)
        ba                  (:,1) double {mustBeVector}
        freq
        ROIs                {mustBeA(ROIs,"ROI")}
        CMfreq              double %nodes x freq
        CMroi               double %nodes x nodes (averaged across frequencies)
        CMavg               double %nodes x nodes x freq, always without nans from non-significant nodes
        filename            string
    end
% condition = 'healthy' or 'PD'
% connectivity_metric = 'coh','icoh','grang'
% connectivitynodes = 186 or 1028
% dataset = 'OMEGA' or 'HCP'
% frequencyrange = 'all' (3:87 or 3:36), [f1 : f2], or 'theta', 'alpha', 'beta low','beta high','beta total','gamma narrow'
% ROIs = 'all', {'M1r','M1l','SMA','preSMA','IFGr','IFGl','mediotemporal left'; 'occipital pole left'} or [X Y Z]
% norm = '_norm' or '' %spatially normalized data
% previous_filename = previous_filename is 

addpath('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\3_TUTTI\rawdataCM')

if conn_nodes == 1018
    symsh = 'Sym';
    symnew = 'Sym';
else
    symsh = '';
    symnew = 'Asym';
end
if smoothening
    smoothening = 'sm20';
else
    smoothening ='';
end
%% determine whether data needs to be re-loaded
filename = [smoothening, 'CMavg', symsh, num2str(conn_nodes),'_',dataset,'_',condition,'_',connectivity_metric,norm,'.mat'];
%% load the correct data

if strcmp(previous_filename, filename) && or(anynan(previous_CMavg) && significant, ~anynan(previous_CMavg) && ~significant)
    load(filename,'frequencies','coorddouble')
    CMavg  = previous_CMavg;
    %CMfreq     = previous_CMfreq;
    %CMroi      = previous_CMroi;
else
    disp(['loading: ' filename])

    
    if isempty(norm)
        load(filename, 'CMavg','frequencies','coorddouble')
        CMavg = CMavg; % nodes x nodes x freq
    else
        load(filename, 'CMavg_norm','frequencies','coorddouble')
        CMavg = CMavg_norm; % nodes x freq
    end

    if isa(CMavg, 'double')
        warning('data is not compressed yet')
    elseif isa(CMavg, 'uint8')
        disp('Decompressing data...')
        CMavg = jvh_convert_uint8_to_doub(CMavg, connectivity_metric);
    else
        error('CMavg has unknown class')
    end

    ROIs  = ROIs.set_coord_target_conn_nodes(conn_nodes);
end

freq  = determine_frequencies(frequencyrange, frequencies);

ROIs  = ROIs.determine_ROI();
%ROIs  = determine_ROIs(ROIs, connectivity_nodes, coorddouble);
%% preset CM

if significant
    load(['ClusterH',symnew, num2str(conn_nodes),'_',dataset,'_',condition,'_',connectivity_metric,norm,'.mat'],'h_cluster')
    temp_ind=find(h_cluster==1);
    temp_vals=CMavg(temp_ind);
    CMavg(h_cluster==0)=NaN; %use "omitnan" in mean function below to still include NaN values as soon as they are accompanied with a non-NaN at point of averaging
end

% need to fix this to know when I need 3D or 2D matrices
if strcmp(norm,'_norm')
    CMfreq=squeeze(CMavg(ROIs.index,:)); %CMfreq is nodes x freq
    CMroi = mean(CMfreq(:,freq.index),2,"omitnan"); %CMroi is here 1 dim
elseif (length(ROIs.index)==1) || ~ROIs.multiple   %(length(ROIs.index)>100)
    CMfreq=squeeze(mean(CMavg(ROIs.index,:,:),1,"omitnan")); %CMfreq is nodes x freq
    CMroi =squeeze(mean(CMavg(:,:,freq.index),3,"omitnan")); %CMroi is nodes x nodes
elseif ROIs.multiple && ~ROIs.restrict
    for i=1:length(ROIs.index)
        CM_inter(:,i,:)=mean(CMavg(ROIs.index,:,:),1,"omitnan"); %CM_inter is nodes x coord_interest x freq
    end
    CMfreq = squeeze(mean(CM_inter,2,"omitnan"));  %CM is nodes x freq
    CMroi = squeeze(mean(CMavg(:,:,freq.index),3,"omitnan")); %CM is nodes x nodes;
elseif ROIs.restricted
        CMfreq=squeeze(CMavg(ROIs.index,ROIs.index,:)); %CM is nodes x nodes x freq => only for small networks
        CMroi=CMfreq; % only for small matrices
else
    error('no ba from CM extracted')
end

ba=mean(CMfreq(:,freq.index),2,"omitnan");
CMfreq = CMfreq(:,freq.index);
CMroi = CMroi(:,ROIs.index);

if significant
    CMavg(temp_ind)=temp_vals;
end

