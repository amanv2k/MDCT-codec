function [out] =  audioDecoder(filename, sparseRowCount, rowCount, colCount, q, Fs)
    fileId = fopen(filename, 'r');
    Raw = fread(fileId, [sparseRowCount, 3], 'int16');
    fclose(fileId);
    n = rowCount;
    nb = colCount;
    A = zeros(n, nb);
    for foo = 1:sparseRowCount
        rowNdx = Raw(foo, 1);
        colNdx = Raw(foo, 2);
        value = Raw(foo, 3);
        A(rowNdx, colNdx) = value;
    end
    M = ones(n, 2*n);       
    for i = 1:n            
        for j = 1:2*n
            M(i, j) = cos((i-1 + 1/2) * (j-1 + 1/2 + n/2) * pi/n);
        end
    end
    N = M';                 
    out = [];
    y1 = A(:, 1);       
    for k=1:nb              
        y1 = A(:, k);       
        y2 = y1*q;          
        w(:, k) = N*y2;        
        if (k>1)            
            w2 = w(n+1:2*n, k-1);
            w3 = w(1:n, k);
            out = [out; (w2 + w3) /2];  
        end
    end
    audiowrite('resconstructed.wav', out, Fs);
    plot(out);  
end