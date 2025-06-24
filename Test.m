function Test(input, numMaps, useLZC, useHoriz, DataHidingKey)
    [~, inputName, ~] = fileparts(char(input));
    logFile = fopen('Results.txt','a');
    channelLabels = ['R', 'G', 'B', 'E'];
    embeddingCapacity = 0;
    fprintf(logFile, '%s ', inputName);
    switch numMaps
        case 4, fprintf(logFile, 'FourMaps ');
        case 1, fprintf(logFile, 'OneMap ');
        case 2, fprintf(logFile, 'TwoMaps ');
    end
    fprintf(logFile, '%s ', ternary(useLZC, 'LZC', 'MSB'));
    rawData = readmatrix(input);
    imgWidth = rawData(1,1);
    imgHeight = rawData(1,2);
    rawData = uint8(rawData(2:end,:));
    binaryLayer = repmat('0', imgHeight, imgWidth, 8);     
    outputImage = zeros(imgHeight, imgWidth, 4, 'uint8');    
    image = zeros(imgHeight, imgWidth, 4, 'uint8');
    embeddedImg = zeros(imgHeight-1, imgWidth-1, 4, 'uint8');
    predictionDiff = zeros(imgHeight-1, imgWidth-1, 4);
    bitLengthMap = zeros(imgHeight-1, imgWidth-1, 4);
    bitCountPerChannel = zeros(1,4);
    embeddingStat = zeros(9, 4);
    secretInput = cell(1,4);
    for idx = 1:size(rawData,1)
        row = ceil(idx / imgWidth);
        col = mod(idx - 1, imgWidth) + 1;
        pixel = rawData(idx, :);
        if max(pixel(1:3)) < 128
            pixel(1:3) = pixel(1:3) * 2;
            pixel(4) = pixel(4) - 1;
        end
        if all(pixel(1:3) == 0)
            pixel(4) = 1;
        end
        image(row, col, :) = pixel;
    end
    %Fig. 7&8(a)
    for ch = 1:4                
        imwrite(uint8(image(:,:,ch)), sprintf('%s_Original_%c.bmp', inputName, channelLabels(ch)));
    end
    reference = extractEdgePixelsToBinary(image);
    for ch = 1:4                
        binBuffer = [];
        for i = 2:imgHeight
            for j = 2:imgWidth
                left = image(i, j-1, ch);
                top = image(i-1, j, ch);
                topLeft = image(i-1, j-1, ch);
                current = image(i, j, ch);
                predicted = MEDPredictor(left, top, topLeft);
                if useLZC
                    [lzVal, bitUsed] = computeLZCDiff(predicted, current);
                    predictionDiff(i-1, j-1, ch) = lzVal;
                    bitLengthMap(i-1, j-1, ch) = bitUsed;
                    binaryLayer(i,j) = abs(int16(predicted) - int16(current));
                    binFragment = dec2bin(binaryLayer(i,j), 8);
                    binBuffer = [binBuffer, binFragment(end - (8 - bitUsed - 1):end)];
                else
                    binCurrent = dec2bin(current, 8);
                    binPred = dec2bin(predicted, 8);
                    binaryLayer(i,j) = current;
                    msbVal = compareBinarySequences(binPred, binCurrent);
                    predictionDiff(i-1, j-1, ch) = msbVal;
                    if msbVal < 9
                        bitUsed = msbVal;
                    else
                        bitUsed = msbVal - 1;
                    end
                    bitLengthMap(i-1, j-1, ch) = bitUsed;
                    binFragment = dec2bin(current, 8);
                    binBuffer = [binBuffer, binFragment(end - (8 - bitUsed - 1):end)];
                end
                bitCountPerChannel(ch) = bitCountPerChannel(ch) + bitUsed;
                embeddingStat(bitUsed, ch) = embeddingStat(bitUsed, ch) + 1;
                if bitUsed >= 4
                    embeddedImg(i-1, j-1, ch) = 255;
                else
                    embeddedImg(i-1, j-1, ch) = 0;
                end
            end
        end
        %Fig. 7&8(b)(c)
        imwrite(uint8(embeddedImg(:,:,ch)), sprintf('%s_EC_Dist_%c_%s.bmp', inputName, channelLabels(ch), ternary(useLZC, 'LZC', 'MSB')));        
        binLen = strlength(binBuffer);
        validLen = floor(binLen / 8) * 8;  
        secretBin = extractBetween(binBuffer, 1, validLen);  
        secretInput{ch} = secretBin;
        %Table 1
        fprintf(logFile, '%c %f ', channelLabels(ch), bitCountPerChannel(ch)/(imgWidth-1)/(imgHeight-1));
        embeddingCapacity = embeddingCapacity + bitCountPerChannel(ch);
    end    
    huffmanTreeStrings = cell(1, 4);
    treeData = 0;
    LabelData = 0;
    if numMaps == 4
        for ch = 1:4
            [encodedBits, treeEncodedSize] = codingGenerate(predictionDiff(:, :, ch), ternary(useLZC, 5, 4));
            treeData = treeData + treeEncodedSize;
            LabelData = LabelData + length(encodedBits) - treeEncodedSize;
            bitChars = char(encodedBits + '0');
            %Table 4
            fprintf(logFile, '%c %f ',channelLabels(ch),(bitCountPerChannel(ch)-length(encodedBits)-length(reference(ch)))/(imgWidth-1)/(imgHeight-1));
            huffmanTreeStrings{ch} = string(bitChars) + reference(ch);
        end
        %Table 2&3
        fprintf(logFile, '%d %d %d \n', treeData, LabelData, treeData+LabelData);
    elseif numMaps == 1 || numMaps == 2
        combinedChannels = predictionDiff(:,:,1);
        if numMaps == 1
            for ch = 2:4
                combinedChannels = cat(1, combinedChannels, predictionDiff(:,:,ch));
            end
        elseif numMaps == 2
            for ch = 2:3
                combinedChannels = cat(1, combinedChannels, predictionDiff(:,:,ch));
            end
        end
        [encodedBits, treeEncodedSize] = codingGenerate(combinedChannels, ternary(useLZC, 5, 4));
        treeData = treeData + treeEncodedSize;
        LabelData = LabelData + length(encodedBits) - treeEncodedSize;
        huffmanTreeStrings{1} = string(encodedBits)';
        if numMaps == 1
            %Table 2&3
            fprintf(logFile, '%d %d %d\n', treeData, LabelData, treeData+LabelData);
        elseif numMaps == 2
	        [encodedBits, treeEncodedSize] = codingGenerate(predictionDiff(:,:,4), ternary(useLZC, 5, 4));
            treeData = treeData + treeEncodedSize;
            LabelData = LabelData + length(encodedBits) - treeEncodedSize;
            %Table 2&3
            fprintf(logFile, '%d %d %d\n', treeData, LabelData, treeData+LabelData);
        end
    end
    if numMaps == 4
        sharedSecret = cell(4,5);
        for ch = 1:4
            bitstream = secretInput{ch};
            segmentCount = floor(strlength(bitstream)/8);
            startIdx = (0:segmentCount-1) * 8 + 1;
            endIdx = startIdx + 7;    
            shareBuffer = cell(1,5);
            for i = 1:5
                shareBuffer{i} = repmat('0', 1, segmentCount * 8);
            end    
            for s = 1:segmentCount
                byte = extractBetween(bitstream, startIdx(s), endIdx(s));
                secretShares = ShamirSecretSharing(3, 5, bin2dec(byte));
                for i = 1:5
                    shareBuffer{i}((s-1)*8+1 : s*8) = dec2bin(secretShares(i,2), 8);
                end
            end    
            for i = 1:5
                sharedSecret{ch,i} = string(shareBuffer{i});
            end
        end
        for shareIndex = 1:5
            rng(DataHidingKey(shareIndex));                   
            for ch = 1:4
                bitIdx = 1;
                bits = char(huffmanTreeStrings{ch});
                for i = 1:imgWidth
                    binaryLayer(1, i, 1:8) = char(bits(bitIdx : bitIdx + 7));
                    bitIdx = bitIdx + 8;
                end
                for i = 2:imgHeight
                    binaryLayer(i, 1, 1:8) = char(bits(bitIdx : bitIdx + 7));
                    bitIdx = bitIdx + 8;
                end
                if useHoriz
                    msbCursor = 1;
                    lsbCursor = 1;                        
                    remainingBits = char(huffmanTreeStrings{ch});
                    remainingBits = remainingBits(bitIdx-1:end);                        
                    huffmanStr = char(remainingBits);
                    secretStr = char(sharedSecret{ch, shareIndex});                
                    for i = 2:imgHeight
                        for j = 2:imgWidth
                            for bit = 1:8
                                if bit <= bitLengthMap(i-1, j-1, ch)
                                    if msbCursor <= strlength(huffmanStr)
                                        binaryLayer(i, j, bit) = huffmanStr(msbCursor); 
                                        msbCursor = msbCursor + 1;
                                    else
                                        binaryLayer(i, j, bit) = char('0' + randi([0 1])); 
                                    end
                                else
                                    if lsbCursor <= strlength(secretStr)
                                        binaryLayer(i, j, bit) = secretStr(lsbCursor);
                                        lsbCursor = lsbCursor + 1;
                                    end
                                end
                            end
                        end
                    end
            
                    for i = 1:imgHeight
                        for j = 1:imgWidth
                            binStr = reshape(binaryLayer(i, j, :), 1, []); 
                            outputImage(i, j, ch) = bin2dec(binStr);
                        end
                    end
                    %Fig. 9
                    imwrite(uint8(outputImage(:,:,ch)), sprintf('%s_Share_%d_%c_%s.bmp', inputName, shareIndex, channelLabels(ch),ternary(useHoriz, 'Horiz', 'Verti')));
                else
                    msbCursor = 1;
                    lsbCursor = 1;                        
                    remainingBits = char(huffmanTreeStrings{ch});
                    remainingBits = remainingBits(bitIdx-1:end);                        
                    huffmanStr = char(remainingBits);
                    secretStr = char(sharedSecret{ch, shareIndex});                
                    for bit = 1:8
                        for i = 2:imgHeight
                            for j = 2:imgWidth
                                if bit <= bitLengthMap(i-1, j-1, ch)
                                    if msbCursor <= strlength(huffmanStr)
                                        binaryLayer(i, j, bit) = huffmanStr(msbCursor);  
                                        msbCursor = msbCursor + 1;
                                    else
                                        binaryLayer(i, j, bit) = char('0' + randi([0 1])); 
                                    end
                                else
                                    if lsbCursor <= strlength(secretStr)
                                        binaryLayer(i, j, bit) = secretStr(lsbCursor);  
                                        lsbCursor = lsbCursor + 1;
                                    end
                                end
                            end
                        end
                    end
                    for i = 1:imgHeight
                        for j = 1:imgWidth
                            binStr = reshape(binaryLayer(i, j, :), 1, []); 
                            outputImage(i, j, ch) = bin2dec(binStr);                                
                        end
                    end
                    %Fig. 9
                    imwrite(uint8(outputImage(:,:,ch)), sprintf('%s_Share_%d_%c_%s.bmp', inputName, shareIndex, channelLabels(ch),ternary(useHoriz, 'Horiz', 'Verti')));
                end
            end
            filename = sprintf('%s_Share_%d_%s.txt', inputName, shareIndex, ternary(useHoriz, 'Horiz', 'Verti'));
            writeRGBEToTxt(outputImage, filename);
        end
    end
    fclose(logFile);
end

function predicted = MEDPredictor(a, b, c)
    if c <= min(a,b)
        predicted = max(a,b);
    elseif c >= max(a,b)
        predicted = min(a,b);
    else
        predicted = int16(a) + int16(b) - int16(c);
    end
end
function [lzVal, count] = computeLZCDiff(pred, actual)
    diff = abs(int16(pred) - int16(actual));
    binStr = dec2bin(diff, 8);
    count = find(binStr ~= '0', 1) - 1;
    if isempty(count), count = 8; end
    if actual >= pred
        lzVal = 8 + count;
    else
        lzVal = 8 - count;
    end
    if count ~= 8
        count = count + 1;
    end
end

function result = ternary(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end

function [encodedBits, treeEncodedSize] = codingGenerate(bytes, symbolSize)
    bytes = bytes(:);  
    uniqueSymbols = unique(bytes);
    freq = histc(bytes, uniqueSymbols);
    symbolCount = length(uniqueSymbols);
    
    nodes = {};
    for i = 1:symbolCount
        nodes{end+1} = struct(...
            'value', uniqueSymbols(i), ...
            'freq', freq(i), ...
            'left', [], ...
            'right', []);
    end

    while length(nodes) > 1
        [~, idx] = sort(cellfun(@(n)n.freq, nodes));
        nodes = nodes(idx);
        left = nodes{1};
        right = nodes{2};
        newNode = struct(...
            'value', [], ...
            'freq', left.freq + right.freq, ...
            'left', left, ...
            'right', right);
        nodes = [{newNode}, nodes(3:end)];
    end
    tree = nodes{1}; 

    codeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    encodedBits = [];  
    stack = {struct('node', tree, 'code', [])};
    
    while ~isempty(stack)
        current = stack{end};
        stack(end) = [];

        if isempty(current.node.left) && isempty(current.node.right)  
            encodedBits(end+1) = true; 
            bits = symbol2bits(current.node.value, symbolSize);
            encodedBits = [encodedBits, bits];
            codeMap(char(current.node.value)) = current.code;
        else
            encodedBits(end+1) = false;  

            stack{end+1} = struct('node', current.node.right, ...
                                  'code', [current.code, false]);
            stack{end+1} = struct('node', current.node.left, ...
                                  'code', [current.code, true]);
        end
    end
    treeEncodedSize = length(encodedBits);
    totalLen = 0;
    for i = 1:length(bytes)
        totalLen = totalLen + length(codeMap(char(bytes(i))));
    end
    
    encodedData = false(1, totalLen);  
    pos = 1;
    for i = 1:length(bytes)
        code = codeMap(char(bytes(i)));
        len = length(code);
        encodedData(pos:pos+len-1) = code;
        pos = pos + len;
    end
    
    encodedBits = [encodedBits, encodedData];
end

function bits = symbol2bits(value, bitLength)
    bits = false(1, bitLength);
    for i = 1:bitLength
        bits(bitLength - i + 1) = bitand(value, 1);
        value = bitshift(value, -1);
    end
end

function binStrArray = extractEdgePixelsToBinary(image)
    [H, W, ~] = size(image);    
    rowPixels = squeeze(image(1, :, :));  
    colPixels = squeeze(image(2:H, 1, :));
    pixels = [rowPixels; colPixels];
    R_bits = "";
    G_bits = "";
    B_bits = "";
    E_bits = "";
    for i = 1:size(pixels, 1)
        R_bits = R_bits + dec2bin(pixels(i, 1), 8);
        G_bits = G_bits + dec2bin(pixels(i, 2), 8);
        B_bits = B_bits + dec2bin(pixels(i, 3), 8);
        E_bits = E_bits + dec2bin(pixels(i, 4), 8);
    end
    binStrArray = [R_bits; G_bits; B_bits; E_bits];
end

function writeRGBEToTxt(outputImage, filename)
    [H, W, C] = size(outputImage);
    fid = fopen(filename, 'w');
    fprintf(fid, '%d %d 0 0\n', W, H);
    for i = 1:H
        for j = 1:W
            pixel = outputImage(i, j, :);
            fprintf(fid, '%d %d %d %d\n', pixel(1), pixel(2), pixel(3), pixel(4));
        end
    end
    fclose(fid);
end
