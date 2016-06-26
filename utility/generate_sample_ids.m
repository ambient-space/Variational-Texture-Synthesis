function sample_ids = generate_sample_ids(sz,dataratio)
    sz=sz(1:2);
    
    N = prod(sz);
    p = randperm(N);
    sample_ids = sort(p(1:ceil(N*dataratio))');
    
end