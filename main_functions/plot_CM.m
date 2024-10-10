function [ba] = plot_CM(cond, conn_metric, conn_nodes, dataset, frequencyrange, ROIs, norm, mode, CMavg, CMfreq, CMroi)

% condition = 'healthy' or 'PD'
% connectivity_metric = 'coh','icoh','grang'
% connectivitynodes = 186 or 1028
% dataset = 'OMEGA' or 'HCP'
% frequencyrange = 'all' (3:87 or 3:36), [f1 : f2], or 'theta', 'alpha', 'beta low','beta high','beta total','gamma narrow'
% ROIs = 'all', {'M1r','M1l','SMA','preSMA','IFGr','IFGl','mediotemporal left'; 'occipital pole left'} or [X Y Z]
% norm = '_norm' or '' %spatially normalized data
% mode = 'whole-brain connectivity', 'spectral plot', 'matrix'
addpath('C:\Users\Jonathan\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT MEG Connectomics\3_TUTTI\rawdataCM')
filename = ['CMavg',num2str(conn_nodes),'_',dataset,'_',cond,'_',conn_metric,norm];
load(filename,'frequencies','coorddouble')
%disp(filename)

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

%need to fix this
colorpanel         = {'black','white'};
%colorpanel         = {'none','black'};
save_indiv_figures = logical(0);

if strcmp(norm,'_norm')
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

% % need to fix this to know when I need 3D or 2D matrices
% if strcmp(norm,'_norm')
%     CM=squeeze(CM(ROI.index,:)); %CM is nodes x freq
% elseif (length(ROI.index)==1) || (length(ROI.index)>100)
%     CM=squeeze(mean(CM(ROI.index,:,:),1)); %CM is nodes x freq
% else
%     CM=squeeze(CM(ROI.index,ROI.index,:)); %CM is nodes x nodes x freq
% end

%% preset low and high color bar
if strcmp(conn_metric,'coh') && strcmp(norm,'_norm')
    con_low=0.5; con_high=1;
elseif strcmp(conn_metric,'coh') && isempty(norm)
    con_low=0; con_high=1;
elseif strcmp(conn_metric,'icoh') && strcmp(norm,'_norm')
    con_low=0.95; con_high=1.05;
elseif strcmp(conn_metric,'icoh') && isempty(norm)
    con_low=0.04; con_high=0.06;
elseif strcmp(conn_metric,'grang') && strcmp(norm,'_norm')
    con_low=0.5; con_high=1.5;
elseif strcmp(conn_metric,'grang') && isempty(norm)
    con_low=-0.1; con_high=0.1;
else
    error('case not found')
end


%% plot figure
if strcmp(mode, 'whole-brain connectivity')
    if strcmp(norm,'_norm') && length(ROI.index) == 1
        error('spatially normalized matrices cannot have whole-brain connectivity with 1 specific region of interest')
    end
	
    ba = plot_whole_brain(freq,CMfreq,coorddouble,ROI, cond, connectivity_metric_label, conn_nodes, mode, dataset, norm_label, con_low, con_high, colorpanel)
elseif strcmp(mode, 'spectral plot')
    assert(length(ROI.index) == 1)
    
    figure
    y=squeeze(mean(CMavg(:,ROI.index,:),1,"omitnan"));
    x=frequencies;
    plot(x,y)
    %ylim([con_low,con_high])
    ylim('auto')
    ba=y;
    title([dataset ' ' cond ' ' num2str(conn_nodes) ' nodes: ' mode ' ' connectivity_metric_label],'color',colorpanel{2})
    subtitle({['ROI: ' ROI.label{1}],norm_label})
    xlabel ('Frequency (Hz)')
    ylabel (connectivity_metric_label)
    set(gcf,'color','white')
    
elseif strcmp(mode, 'matrix')
 
    if strcmp(norm,'_norm')
        error('spatially normalized matrices have no graph'); %CM is nodes x freq
    end
    
    if length(freq.frequencybandnames) > 1
        figure
        set(gca,'Position',[0.1 0.1 0.9 0.9],'units','normalized','outerposition',[0 0 1 1])
        sgtitle({[dataset ' ' cond ' ' num2str(conn_nodes) ' nodes: ' mode ' ' connectivity_metric_label],['ROI: ' ROI.label{1}],norm_label},'color',colorpanel{2})  
        for i=1:length(freq.frequencybandnames)
            subplot(1,length(freq.frequencybandnames),i)
            figone(8,50)
            CMfreq=squeeze(CMavg(:,:,freq.frequency_indices{i}));
            CMfreq=mean(CMfreq,3,"omitnan");
            imagesc(CMfreq)
            title(freq.frequencybandnames{i},'color',colorpanel{2}, 'interprete','none','fontsize',10,'Margin',1 )

            axis equal
            axis off
            set(gcf,'color',colorpanel{1})
            set(gca,'color',colorpanel{2},'ycolor',colorpanel{2},'xcolor',colorpanel{2})

            hold off
        end

        hold on
        c = colorbar('eastoutside');
        c.Color =colorpanel{2};

        ba = CMfreq;
    else
        figure
        CMfreq=squeeze(CMavg(:,:,freq.frequency_indices{1}));
        CMfreq=mean(CMfreq,3,'omitnan');
        imagesc(CMfreq,[con_low con_high])
        %xticklabels(ROI.label)
        %yticklabels(ROI.label)
        colorbar
        title({[dataset ' ' cond ' ' num2str(conn_nodes) ' nodes: ' mode ' ' connectivity_metric_label]},'color',colorpanel{2})  
        ba = CMfreq;
        set(gcf,'color',colorpanel{1})
        set(gca,'color',colorpanel{2},'ycolor',colorpanel{2},'xcolor',colorpanel{2})
    end
end

if save_indiv_figures; save_indiv_figures_func(); end


end

function ba = plot_whole_brain(freq,CM,coorddouble,ROI, condition, connectivity_metric_label, connectivity_nodes, mode, dataset, norm_label, con_low, con_high, colorpanel)
        
    ctx=export(gifti('BrainMesh_ICBM152_smoothed.gii'));
    
    figure
    set(gca,'Position',[0.1 0.1 0.9 0.9],'units','normalized','outerposition',[0 0 1 1])
    for a=1:size(coorddouble,1)
        [mind(1,a),i(1,a)] = min(wjn_distance(ctx.vertices,[coorddouble(a,1:3)]));
        nmni(a,:) = ctx.vertices(i(a),:);
    end

    if length(freq.frequencybandnames) == 6
        subplot_total = 6;
    elseif length(freq.frequencybandnames) == 5
        subplot_total = 5;
    else
        subplot_total = 4;
    end

    for i=1:length(freq.frequencybandnames)
        ba = mean(CM(:,freq.frequency_indices{i}),2);
        subplot(3,subplot_total,[i,i+subplot_total])
        p=wjn_plot_surface(ctx);
        figone(15,50)
        view(0,90)
        camlight 
        material dull
        hold on
        cm = colormap('jet');
        [s,v,ov]=jvh_plot_colored_spheres(nmni,ba,3,cm,con_low,con_high, ROI.index);
        title(freq.frequencybandnames{i},'color',colorpanel{2}, 'interprete','none','fontsize',10,'Margin',1)
        hold off
    end
    if i == 1
        i=2:3;
        subplot(3,subplot_total,[i,i+subplot_total])
        p=wjn_plot_surface(ctx);
        %zoom(1.05)
        view(-90,-15)
        camlight 
        material dull
        hold on
        cm = colormap('jet');
        [s,v,ov]=jvh_plot_colored_spheres(nmni,ba,3,cm,con_low,con_high, ROI.index);
        title(freq.frequencybandnames{1},'color',colorpanel{2}, 'interprete','none','fontsize',10,'Margin',1)
        hold off
        i=4;
        subplot(3,subplot_total,[i,i+subplot_total])
        p=wjn_plot_surface(ctx);
        zoom(0.8)
        %figone(40,40)
        view(0,0)
        camlight 
        material dull
        hold on
        cm = colormap('jet');
        [s,v,ov]=jvh_plot_colored_spheres(nmni,ba,3,cm,con_low,con_high, ROI.index);
        title(freq.frequencybandnames{i-3},'color',colorpanel{2}, 'interprete','none','fontsize',10,'Margin',1)
        hold off
    end

    hold on
    subplot(3,subplot_total,[3*subplot_total-subplot_total+1:3*subplot_total])
    c = colorbar('North');
    c.Color =colorpanel{2};
    limits = [con_low,con_high];
    set(gca,'clim',limits([1,end]))
    axis('off')
    sgtitle({[dataset ' ' condition ' ' num2str(connectivity_nodes) ' nodes: ' mode ' ' connectivity_metric_label],['ROI: ' ROI.label{1}],norm_label},'Color',colorpanel{2})
    hold off
    set(gcf,'color',colorpanel{1})
end

% function ROI  = determine_ROIs(ROIs, connectivitynodes, coorddouble)
%     % ROIs = 'all', {'M1r','M1l','SMA','preSMA','IFGr','IFGl','mediotemporal left'; 'occipital pole left'} or [X Y Z]
%     if isnumeric(ROIs)
%         if isequal(size(ROIs),[1 3])
%             ROI.MNI        = ROIs;
%             ROI.label      = {'custom'};
%             [ROI.target, ROI.index] = jvh_mindistancecoords(ROI.MNI,coorddouble);
%         elseif isequal(size(ROIs,2),3)
%             ROI.MNI        = ROIs;
%             error('test motor stopping network')
%             for el=1:size(ROIs,1)
%                 [ROI.target(el,:), ROI.index(el)] = jvh_mindistancecoords(ROI.MNI(el,:),coorddouble);
%                 ROI.label      = {'SMAr','SMAl','preSMAr','preSMAl','IFGr','IFGl'}
%             end
%         else
%             error('no customized MNI could be retrieved')
%         end
%     elseif any(ismember({'M1r','M1l','SMAr','SMAl','preSMAr','preSMAl','IFGr','IFGl', 'mediotemporal left', 'occipital pole left'},ROIs))     
%         names       = {'M1r','M1l','SMAr','SMAl','preSMAr','preSMAl','IFGr','IFGl','mediotemporal left', 'occipital pole left'};
%         mni         = [37   -25    64; -37   -25    62;  2   -10    59; -2   -10    59; 2    30    48;  -2    30    48;  42    26    14;  -42    26    14; -61 -14 -13; -15 -99 12];
%         ROI.MNI                 = mni(ismember(names, ROIs),:);
%         ROI.label               = {ROIs};
%         [ROI.target, ROI.index] = jvh_mindistancecoords(ROI.MNI,coorddouble);
%     elseif strcmp(ROIs, 'all')
%         ROI.label = {'all'};
%         ROI.index = [1:connectivitynodes];
%     elseif isfile(ROIs)
%         mask_nii = ea_load_nii(ROIs)
%         %calculate the centroid point of the mask and rearrange dimensions
%         measurement= regionprops(mask_nii.img, 'Centroid');
%         measurement.Centroid = [measurement.Centroid(2),measurement.Centroid(1),measurement.Centroid(3)]; %put it in x y z
% 
%         %convert coordinate to MNI space for the ROI
%         ROI.MNI = wjn_cor2mni(measurement.Centroid, mask_nii.mat);
%         ROI.label = {'mask'};
%         %calculate distance between centroid of mask and closest cortical node
%         [ROI.target, ROI.index] = jvh_mindistancecoords(ROI.MNI,coorddouble);
%     else
%         error('ROI region not found')
%     end
% 
% end

% function [freq]=determine_frequencies(frequencyrange, frequencies)
%     % frequencyrange = 'all' (3:87 or 3:36), [f1 : f2], or 'theta', 'alpha', 'beta_low','beta_high','beta_total','gamma_narrow'
% if nargin < 2
%     load('standard_frequencies.mat','frequencies')
% end
% if isnumeric(frequencyrange)
%     assert(frequencyrange(1)>= min(frequencies))
%     assert(frequencyrange(2)<= max(frequencies))
%     freq.frequency_indices{1} = frequencyrange - 3;
%     freq.frequencybandnames={'custom'};
%     freq.frequencyband{1} = frequencyrange;
% elseif ischar(frequencyrange)
%     freq.frequencybandnames={'theta', 'alpha', 'beta low','beta high','beta total'};
%     freq.frequencyband{1} = 4:7;
%     freq.frequencyband{2} = 7:12;
%     freq.frequencyband{3} = 13:21;
%     freq.frequencyband{4} = 21:35;
%     freq.frequencyband{5} = 13:35;
%     if max(frequencies) == 90
%         freq.frequencybandnames{6} = 'gamma narrow';
%         freq.frequencyband{6} = 60:90;        
%     end
%     for i = 1:length(freq.frequencyband)
%         freq.frequency_indices{i} = freq.frequencyband{i}-3;
%     end
%     if strcmp(frequencyrange,'all')
%         return
%     else
%         ind=ismember(freq.frequencybandnames,frequencyrange);
%         freq.frequencybandnames = {frequencyrange};
%         temp   = freq.frequencyband{ind};
%         freq= rmfield(freq, 'frequencyband');
%         freq.frequencyband{1}   = temp;
%         temp   = freq.frequency_indices{ind};
%         freq= rmfield(freq, 'frequency_indices');
%         freq.frequency_indices{1}   = temp;
%     end
% 
% else
%     error("frequencyrange undefined")
% end
% 
% 
% end

function save_indiv_figures_func()
    disp('printing figure now')
    print(gcf,[datestr(now,'yyyymmdd_HHMMSS') '.png'],'-dpng','-r2000'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_2_hemisphere()
    wjn_spherical_roi('SMA.nii',[-2 -10 59],3*2.5)
    spm_surf('SMA.nii',2)
    wjn_spherical_roi('preSMA.nii',[2 30 48],3*2.5)
    spm_surf('preSMA.nii',2)
    wjn_spherical_roi('M1r.nii',[37 -25 64],3*2.5)
    spm_surf('M1r.nii',2)
end

function calculate_difference_PC_HC()
    connectivitynodes = 186;
    
    load('CMavg186_OMEGA_healthy_coh.mat')
    ROI = determine_ROIs({'SMA'}, connectivitynodes, coorddouble);
    freq = determine_frequencies('beta total', frequencies);
    CM = CMavg;
    CM = squeeze(mean(CM(ROI.index,:,:),1));
    ba = mean(CM(:,freq.frequency_indices{1}),2);
    HC = ba;
    load('CMavg186_OMEGA_PD_coh.mat')
    CM = CMavg;
    CM = squeeze(mean(CM(ROI.index,:,:),1));
    ba = mean(CM(:,freq.frequency_indices{1}),2);
    PD = ba;
    
    conn = PD - HC;
    wjn_heatmap('PD-HC_OMEGA_coh_beta_SMA.nii',coorddouble,conn,'HarvardOxford-cort-maxprob-thr0-1mm.nii')
    
  
    load('CMavg186_OMEGA_healthy_coh.mat')
    ROI = determine_ROIs({'preSMA'}, connectivitynodes, coorddouble);
    freq = determine_frequencies('theta', frequencies);
    CM = CMavg;
    CM = squeeze(mean(CM(ROI.index,:,:),1));
    ba = mean(CM(:,freq.frequency_indices{1}),2);
    HC = ba;
    load('CMavg186_OMEGA_PD_coh.mat')
    CM = CMavg;
    CM = squeeze(mean(CM(ROI.index,:,:),1));
    ba = mean(CM(:,freq.frequency_indices{1}),2);
    PD = ba;
    
    conn = PD - HC;
    wjn_heatmap('PD-HC_OMEGA_coh_theta_preSMA.nii',coorddouble,conn,'HarvardOxford-cort-maxprob-thr0-1mm.nii')
    
    load('CMavg186_OMEGA_healthy_coh.mat')
    ROI = determine_ROIs({'M1r'}, connectivitynodes, coorddouble);
    freq = determine_frequencies('beta total', frequencies);
    CM = CMavg;
    CM = squeeze(mean(CM(ROI.index,:,:),1));
    ba = mean(CM(:,freq.frequency_indices{1}),2);
    HC = ba;
    load('CMavg186_OMEGA_PD_coh.mat')
    CM = CMavg;
    CM = squeeze(mean(CM(ROI.index,:,:),1));
    ba = mean(CM(:,freq.frequency_indices{1}),2);
    PD = ba;
    
    conn = PD - HC;
    wjn_heatmap('PD-HC_OMEGA_coh_beta_M1r.nii',coorddouble,conn,'HarvardOxford-cort-maxprob-thr0-1mm.nii')

end