function sample_ids = generate_sample_ids(sz,blocksize,dataratio)
    sz=sz(1:2);
%     samplematrix is initialized withs 1s at sample locations
    samplematrix = false(sz);
    SampleIndex = randperm(prod(sz));
    samplematrix(SampleIndex(1: fix(dataratio*numel(samplematrix)))) = true;
    %patches are indexed at topleft of patch; pels within the
    %last patch at the right/bottom edge of image can't be sampled
    samplematrix(end-blocksize+2:end,:)=false;
    samplematrix(:,end-blocksize+2:end)=false;
    sample_ids = find(samplematrix);
end