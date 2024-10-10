%% This function converts an cortical brain atlas to ConnectoMEG compatible mesh
%
% atlaslabelstext: is csv or text file with in the first column the numbers of the parcellations, in the second column the label of the parcellation
% atlasnifti: the nifti file (.nii or .nii.gz) with the atlas, with corresponding numbers of the parcellation
% atlasname: short alphanumeric name for saving purpose (exclude white space or special characters)
% exclude_areas: a vector with numbers of areas to exclude (such as subcortical areas)
%
% Code to generate mesh2labels:
% add_atlas_to_mesh2labels('AAL2.txt','AAL2.nii.gz','AAL2',[83,84,95:120])
% add_atlas_to_mesh2labels('Brainnectome.txt','Brainnectome.nii.gz','Brainnectome',[211:246])
% add_atlas_to_mesh2labels('Yeo2011_7Networks_MNI152 (Yeo 2011).txt','Yeo2011_7Networks_MNI152 (Yeo 2011).nii','Yeo7')
% add_atlas_to_mesh2labels('Mindboggle 101 - Desikan protocol (Klein 2012).txt','Mindboggle 101 - Desikan protocol (Klein 2012).nii','Mindboggle',[1:1001])

function mesh2labels = add_atlas_to_mesh2labels(atlaslabelstext, atlasnifti, atlasname, exclude_areas)
% addpath(genpath('C:\Users\Jonathan\Documents\MATLAB\add_on_Matlab\lead'))
% addpath('C:\Users\Jonathan\Documents\MATLAB\add_on_Matlab\spm12')
% 
arguments
    atlaslabelstext {mustBeFile}
    atlasnifti {mustBeFile}
    atlasname {mustBeTextScalar}
    exclude_areas {mustBeVector} = [0]
end
    [~, ~, fExt] = fileparts(atlaslabelstext);
    assert(ismember(fExt,{'.txt','.csv'}))
    [~, ~, fExt] = fileparts(atlasnifti);
    assert(ismember(fExt,{'.nii','.gz'}))
    assert(regexp(atlasname, '^[A-Za-z0-9]+$'))
    warning('lead dbs is needed for this function')

    T = readtable(atlaslabelstext, 'Delimiter',' ');
    T = rmmissing(T,2);
    uitable('Data',table2cell(T(:,:)),'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    assert(width(T)>1,'the table needs to have space separated columns: column 1 is parcel number, column 2 the parcel label')
    MNI_template_labels = ea_load_nii(atlasnifti);
    %MNI_template_info = niftiinfo(atlasnifti);
    
    if any(exclude_areas)
        for x=exclude_areas
               MNI_template_labels.img(MNI_template_labels.img == x) = 0;
               T(T.(1)==x,:)=[];
        end
    end
    
    [cor_x,cor_y,cor_z]=ind2sub(size(MNI_template_labels.img),find(MNI_template_labels.img>0));  %check for labels with voxel value more than 0
    labels_values=MNI_template_labels.img(sub2ind(size(MNI_template_labels.img),cor_x,cor_y,cor_z));
    labels_mni=cor2mni([cor_x,cor_y,cor_z],MNI_template_labels.mat);
    path_to_save = which('mesh2labels.mat');
    load('mesh2labels.mat')
    coorddoublename={'coorddouble186.mat','coorddouble1018.mat','coorddouble1028.mat'};
    coordsymmetry={'asym186','sym1018','asym1028'};
    number_of_nodes=[186,1018,1028];
    for c=1:3
        load(coorddoublename{c}, 'coorddouble')
        [closest_coord, index, distance]=jvh_mindistancecoords(coorddouble,labels_mni);
        
        coorddouble_labelnr=labels_values(index);
        
        myLabel={};
        for i=1:number_of_nodes(c)
            myLabel{i}=T.(2){T.(1)==coorddouble_labelnr(i)};
        end
        myLabel=replace(myLabel,'_',' ');
        
        mesh2labels.([coordsymmetry{c} '_to_' atlasname]).coorddouble_labelnr = coorddouble_labelnr;
        mesh2labels.([coordsymmetry{c} '_to_' atlasname]).labelnr_unique = T.(1);
        mesh2labels.([coordsymmetry{c} '_to_' atlasname]).myLabel = myLabel';
        mesh2labels.([coordsymmetry{c} '_to_' atlasname]).myLabel_unique = T.(2);
        
        mesh2labels.Items(end+1) = {atlasname};
        save(path_to_save,"mesh2labels")
    end

end

function mni = cor2mni(cor, T)
    % function mni = cor2mni(cor, T)
    % convert matrix coordinate to mni coordinate
    %
    % cor: an Nx3 matrix
    % T: rotation matrix
    % mni is the returned coordinate in mni space
    %
    % xu cui
    % 2004-8-18
    % last revised: 2005-04-30
    
    % if nargin == 1
    %     T = ...
    %         [-4     0     0    84;...
    %          0     4     0  -116;...
    %          0     0     4   -56;...
    %          0     0     0     1];
    % end
    
    cor = round(cor);
    mni = T*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
    mni = mni';
    mni(:,4) = [];
end