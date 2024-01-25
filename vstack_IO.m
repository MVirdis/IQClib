function M = vstack_IO(A,B,C,D)
%VSTACK_IO returns [I,0;A,B;C,D] given A, B, C, D

M = [eye(size(A,2)), zeros(size(A,2), size(B,2)); A,B;C,D];

end
