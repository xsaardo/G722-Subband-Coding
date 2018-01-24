function out = addTransmissionNoise(IL,BER)
% For applying transmission noise

% Convert codeword to binary representation
binRep = dec2bin(IL,6);
logBinRep = zeros(size(binRep));
logBinRep(binRep == '1') = 1;

% Randomly flip bits according to specified BER
flipchance = (rand(size(binRep)) < BER);
corruptedBin = xor(logBinRep,flipchance);

% Convert from binary representation to codeword
corruptedCell = cell(size(corruptedBin));
corruptedCell(corruptedBin == 1) = {'1'};
corruptedCell(corruptedBin == 0) = {'0'};
out = zeros(size(corruptedCell,1),1);
for i = 1:size(corruptedCell,1)
    out(i) = bin2dec(strjoin(corruptedCell(i,:)));
end

end