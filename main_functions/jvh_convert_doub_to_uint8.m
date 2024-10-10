function int=jvh_convert_doub_to_uint8(doub,conn_metric)
    arguments
        doub {mustBeA(doub,"double")}  
        conn_metric {mustBeMember(conn_metric,['coh','icoh','granger'])}
    end
    switch conn_metric
        case {'coh','icoh'}
            int=uint8(rescale(abs(doub), 0, 255, 'InputMin', 0, 'InputMax', 1));
        case {'grang','rcgranger'}
            int=uint8(rescale(doub, 0, 255, 'InputMin', -25, 'InputMax', 25));
        otherwise
            error('no convertion possible')
    end
end