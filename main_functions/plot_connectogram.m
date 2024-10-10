 %% test connectogram 2 based on percentile
 function plot_connectogram(CM, mesh2labels, percentile)
     arguments
         CM double
         mesh2labels {mustBeA(mesh2labels,"struct")}
         percentile {mustBeInRange(percentile, 0,100, 'exclusive')} = 99;
     end
    % Create custom colormap
    figure                 % Creates a figure
    set(gca,'FontSize',4) % Creates an axes and sets its FontSize to 18
    
    addpath('C:\Users\Jonathan\Documents/CODE/circularGraph-master/circularGraph-master/')
    %load('mesh2labels.mat','mesh2labels') %from jvh toolbox

    index_sorted=mesh2labels.index_sorted;
    myLabels_sorted=mesh2labels.myLabels_sorted;
    myColorMap=mesh2labels.myColorMap;
    %load('CMavg186_OMEGA_PD_coh.mat')

    %CMavg=mean(CMavg(:,:,21-3:35-3),3);
    
    %CMavg=squeeze(CMavg);
    CM(CM==1)=0;
    CM(isnan(CM))=0;
    cutoff = prctile(CM(:),percentile);
    CM(CM>=cutoff)=1;
    CM(CM<cutoff)=0;
    
    %x2= CMavg;
    if sum(CM,'all')>1000
        error('the percentile is too low')
    end
    x = CM; %= abs(x2 - x1);

    x_sorted = x(index_sorted,:);
    CIR=circularGraph(x_sorted,'Colormap',myColorMap(index_sorted,:),'Label',myLabels_sorted);
    M=num2cell(myColorMap,2);
    [CIR.Node.Color] = deal(M{:});
    camroll(-90)
    H=findobj(gca,'Type','text');
    for i=1:186
        H(i).Rotation = H(i).Rotation-90;
    end
    H=findobj(gca,'Type','text');
    for i=1:186
        if ismember(H(i).String,{'left limbic','left frontoparietal','left parietal','left occipital'})
            H(i).HorizontalAlignment = 'right';
            H(i).Rotation = H(i).Rotation-180;
        elseif ismember(H(i).String,{'right frontal','right temporal'})
            H(i).HorizontalAlignment = 'left';
            H(i).Rotation = H(i).Rotation-180;
        end
    end
end



