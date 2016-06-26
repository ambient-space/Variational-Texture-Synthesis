function read_params(ob,params)
    
    %transfer any params into an object
    
    fields = fieldnames(params);
    
    for k = 1:numel(fields)
        
        fn = fields{k};
        if numel(findprop(ob,fn))
            eval(['ob.',fn,' = params.',fn,';']);
        else
            error(['Parameter field ''',fn,''' is not a property of the class ''',class(ob),'''']);
        end
    end
end