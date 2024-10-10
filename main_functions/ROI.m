classdef ROI
    properties
        index           {mustBeNumericOrLogical, mustBeVector} = false;
        coord_closest   (:,3)double {mustBeNumericOrLogical}
        coord_target    (:,3)double {mustBeNumeric}
        coord_interest  (:,3)double {mustBeNumericOrLogical,mustBeInRange(coord_interest,-120,120)}
        label           {mustBeVector,mustBeA(label,"cell")} = {'custom'};
        mask            = '';
        conn_nodes      {mustBeMember(conn_nodes,[186,1018,1028])} = 186;
        multiple        {mustBeNumericOrLogical} = false;
        restrict        {mustBeNumericOrLogical} = false;
    end
    methods
        function obj = ROI()
            obj=set_coord_target_conn_nodes(obj,obj.conn_nodes);
        end
        function obj = set_coord_interest(obj,val)
            if nargin == 2
                obj.coord_interest = val;
                obj = obj.set_label('custom');
                if size(val,1)>1
                    obj.multiple = true;
                else
                    obj.multiple = false;
                end
            end
        end
        function obj = add_coord_interest(obj,val)
            if nargin == 2
                obj.coord_interest = [obj.coord_interest;val];
                if length(obj.label) ~= size(obj.coord_interest,1)
                    obj = obj.set_label('custom');
                end
                obj.multiple = true;
            end
        end
        function obj = set_label(obj,val)
            obj.label = {val};
            if any(ismember({'M1r','M1l','SMAr','SMAl','preSMAr','preSMAl','IFGr','IFGl', 'mediotemporal right','mediotemporal left', 'occipital pole right','occipital pole left'},val))
                names       = {'M1r','M1l','SMAr','SMAl','preSMAr','preSMAl','IFGr','IFGl', 'mediotemporal right','mediotemporal left', 'occipital pole right','occipital pole left'};
                mni         = [37   -25    64; -37   -25    62;  2   -10    59; -2   -10    59; 2    30    48;  -2    30    48;  42    26    14;  -42    26    14; 61 -14 -13; -61 -14 -13; 15 -99 12; -15 -99 12];
                obj.coord_interest   = mni(ismember(names, val),:);
            elseif strcmp(val, 'all')
                obj.index = 1:obj.conn_nodes;
            end
        end
        function obj = add_label(obj,val)
            obj.label(end+1)= {val};
            if any(ismember({'M1r','M1l','SMAr','SMAl','preSMAr','preSMAl','IFGr','IFGl', 'mediotemporal right','mediotemporal left', 'occipital pole right','occipital pole left'},val))
                names       = {'M1r','M1l','SMAr','SMAl','preSMAr','preSMAl','IFGr','IFGl', 'mediotemporal right','mediotemporal left', 'occipital pole right','occipital pole left'};
                mni         = [37   -25    64; -37   -25    62;  2   -10    59; -2   -10    59; 2    30    48;  -2    30    48;  42    26    14;  -42    26    14; 61 -14 -13; -61 -14 -13; 15 -99 12; -15 -99 12];
                obj =obj.add_coord_interest(mni(ismember(names, val),:));
            else
                warning('coordinates not updated')
            end
        end
        function obj = switch_restriction(obj,OnOff)
            if strcmp(OnOff, 'On')
                obj.restrict = 1;
            elseif strcmp(OnOff, 'Off')
                obj.restrict = 0;
            end
        end
        function obj = set_mask(obj,val, path)
            arguments
                obj
                val
                path='';
            end
            if path
                cwd = cd;
                cd(path)
            end
            obj.mask = val;
            if exist("ea_load_nii", "file")
                nii = ea_load_nii(obj.mask);
                [cor_x,cor_y,cor_z]=ind2sub(size(nii.img), find(nii.img>0));
                obj =  obj.set_coord_interest(ROIcor2mni(cor_x,cor_y,cor_z,nii.mat));
            else
                error('add lead dbs to the path')
            end
            if path
                cd(cwd)
            end
        end
        function obj = set_index(obj,val)
            obj.index = val;
            obj.coord_closest = [false false false];
            obj.coord_interest = [false false false];
            if length(val)>1
                obj.multiple=true;
            else
                obj.multiple=false;
            end
            obj=obj.set_label('custom');
        end
        function obj = set_coord_target_conn_nodes(obj,val)
            load(['coorddouble' num2str(val)],'coorddouble')
            obj.coord_target = coorddouble;
            obj.conn_nodes = val;
        end
        function obj = determine_ROI(obj)
            if strcmp(obj.label, 'all')
                obj.index = 1:obj.conn_nodes;
            elseif any(obj.coord_interest)
                [obj.coord_closest, obj.index] = jvh_mindistancecoords(obj.coord_interest, obj.coord_target);
                obj.coord_closest=unique(obj.coord_closest,'rows');
                obj.index=unique(obj.index);
            elseif ~any(obj.coord_interest)
                obj.coord_interest = obj.coord_target(obj.index,:);
            else
                error('ROI mech not yet implemented')
            end
        end
    end
end

% function mustBeValidMask(filename)
%     if ischar(filename)
%         if any(endsWith(filename,[".nii.gz",".nii"]))  
%             mustBeFile(filename)
%         elseif isempty(filename)
%         else
%             error('extension of file is incorrect')
%         end
%     else
%         error('filename must be a string')
%     end
% end

function mni=ROIcor2mni(cor_x, cor_y,cor_z, T)
    mni = T*[cor_x cor_y cor_z ones(size(cor_x,1),1)]';
    mni = mni';
    mni(:,4) = [];
end
