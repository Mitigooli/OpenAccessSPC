function area = overlap( A, B )

area = 0;
if ~( A(3) < B(1) || B(3) < A(1) || A(4) < B(2) || B(4) < A(2) )
	area = 1;
end
