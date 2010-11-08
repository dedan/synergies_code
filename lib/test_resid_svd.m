function v=test_resid_svd( matin )

% simple residual test


SS = sum(sum(matin));
N2 = 20;
v  = NaN(1,N2);
for i=1:N2
   [w,h,score]=my_svd(matin',i);
   v(i) = score;
end
v = (v/SS) * 100;
