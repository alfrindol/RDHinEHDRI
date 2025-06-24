function index = compareBinarySequences(seq1, seq2)    
    if length(seq1) ~= length(seq2)
        error('Sequences must have the same length.');
    end

    for i = 1:length(seq1)
        if seq1(i) ~= seq2(i)
            index = i;
            return;
        end
    end
        
    index = 9;
end