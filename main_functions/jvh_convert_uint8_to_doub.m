function doub=jvh_convert_uint8_to_doub(int,conn_metric)
    arguments
        int {mustBeA(int,"uint8")}  
        conn_metric {mustBeMember(conn_metric,['coh','icoh','granger'])}
    end
    switch conn_metric
        case {'coh','icoh'}
            doub=rescale(double(int), 0, 1, 'InputMin', 0, 'InputMax', 255);
        case {'grang','rcgranger'}
            doub=rescale(double(int), -25, 25, 'InputMin', 0, 'InputMax', 255);
        otherwise
            error('no convertion possible')
    end
end