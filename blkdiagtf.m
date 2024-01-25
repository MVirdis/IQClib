function T = blkdiagtf(varargin)
%BLKDIAGTF diagonal concatenation of transfer matrices

T = varargin{1};
[rows_,cols_] = size(T);
for i=2:nargin
    [rows,cols] = size(varargin{i});
    T = [T, zeros(rows_,cols);
         zeros(rows,cols_), varargin{i}];
    rows_ = rows_ + rows;
    cols_ = cols_ + cols;
end

end
