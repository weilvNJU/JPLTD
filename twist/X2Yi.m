function Y = X2Yi(X, i)

Y = Trans_Faces(shiftdim(X, i-1));  %N*N*V ��V*N*N��N*V*N 

% aaaa = shiftdim(T_tensor,2);   �� T_tensor: a*b*c,�� aaaΪ c*a*b
% B = shiftdim(A,n) ������ A ��ά���ƶ� n ��λ�á�
% �� n Ϊ������ʱ��shiftdim �����ƶ�ά�ȣ�
% �� n Ϊ������ʱ�������ƶ�ά�ȡ�
% ���磬��� A �� 2��3��4 ���飬�� shiftdim(A,2) ���� 4��2��3 ���顣