%Aman Verma 19BEC1284
%Ayush Singh 19BEC1032
%Shreyansh Kumar 19BEC1246
%Parth Desai 19BEC1351
%Topic- Implementation of encoder and decoder on an audio file
%Purpose- To provide secure communication between two persons
function [Fs, q, x, rowCount, colCount, sparseRowCount, rowVector, colVector, v] =  audioEncoder(infilename, outfilename)
    [x, Fs] = audioread(infilename);
    display(Fs)
    plot(x)
    n = 32;             % length of Frame
    nb = floor(length(x)/n - 1);           % number of Frame, > 1
    rowCount = n
    colCount = nb
    b = 4; L = 5;       % Quantization information
    q = 2*L/(2^b - 1)   % b bits on the interval [-L, L]
    for i = 1:n
        for j = 1:2*n
            M(i, j) = cos((i-1 + 1/2) * (j-1 + 1/2 + n/2) * pi/n);
        end
    end
    M = sqrt(2/n) * M; 
    W = ones (n, nb);
    for k=1:nb          % loop over each window
        x0 = x(1+(k-1)*n : 2*n+(k-1)*n);
        y0 = M*x0;
        y1 = round(y0/q); % transform components quantized
        W(:, k) = y1;
    end
    
    [rowVector, colVector, v] = find(W);
    [sparseRowCount, ~] = size(rowVector)
    fileId = fopen(outfilename, 'w');   % open the output file with write access
    fwrite(fileId, [rowVector colVector v], 'int16');         % write the matrix to the file
    fclose(fileId);                     % close the output file
end