function [Fs, q, x, rowCount, colCount, sparseRowCount, rowVector, colVector, v] =  main(infilename, outfilename)
    [x, Fs] = audioread(infilename);
    coef=mdct(x,kbdwin());
%     n = 32;             % length of Frame
%     nb = floor(length(x) / n - 1);           % number of Frame, > 1
%     rowCount = n;
%     colCount = nb;
%     b = 4; L = 5;       % Quantization information
%     q = 2*L/(2^b - 1)  % b bits on the interval [-L, L]
%     for i = 1:n
%         for j = 1:2*n
%             M(i, j) = cos((i-1 + 1/2) * (j-1 + 1/2 + n/2) * pi/n);
%         end
%     end
%     M = sqrt(2/n) * M; 
%     W = ones (n, nb);
%     for k=1:nb          % loop over each window
%         x0 = x(1+(k-1)*n : 2*n+(k-1)*n);
%         y0 = M*x0;
%         y1 = y0/q; % transform components quantized
%         W(:, k) = y1;
%     end
    
    
%     [rowVector colVector v] = find(W);
%     [sparseRowCount, ~] = size(rowVector)
    fileId = fopen(outfilename, 'w');   % open the output file with write access
    fwrite(fileId, coef);         % write the matrix to the file    
    fclose(fileId);                     % close the output file
    %decoding below
    %dlmwrite('sparsetestDLM', [row, col, v], 'delimiter', '\t');

    %A = dlmread(filename);  % read the compressed file into a matrix, A
    %[n, nb] = size(A);      % use the dimensions of the matrix to determine n and nb
    % Raw will be the [row col value] matrix of a sparse matrix
    fileId = fopen(outfilename, 'r');
    Raw=fread(fileId)
    out=imdct(Raw,kbdwin(512))
%     Raw = fread(fileId, [sparseRowCount, 3], 'int16');
%     fclose(fileId);
%     n = rowCount;
%     nb = colCount;
%     A = zeros(n, nb);
%     for foo = 1:sparseRowCount
%         rowNdx = Raw(foo, 1);
%         colNdx = Raw(foo, 2);
%         value = Raw(foo, 3);
%         A(rowNdx, colNdx) = value;
%     end
%     M = ones(n, 2*n);       % initialize our DCT matrix
%     for i = 1:n             % Construct our DCT matrix
%         for j = 1:2*n
%             M(i, j) = cos((i-1 + 1/2) * (j-1 + 1/2 + n/2) * pi/n);
%         end
%     end
%     N = M';                 % since M is orthogonal, its inverse is its transpose
%     out = [];
%     y1 = A(:, 1);
%     for k=1:nb              % loop over each window
%         y1 = A(:, k);       % grab the kth column of A
%         y2 = y1*q;          % dequantize it
%         w(:, k) = N*y2;     % invert the MDCT    
%         if (k>1)            % if we're not on the first step
%             w2 = w(n+1:2*n, k-1);
%             w3 = w(1:n, k);
%             out = [out; (w2 + w3) /2];  % collect the reconstructed signal
%         end
%     end
    audiowrite('resconstructed.wav', out, Fs);
    %sound(out, Fs);  
    plot(out)
end