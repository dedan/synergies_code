function v=test_resid_nmf( matin , Nitr)

SS = sum(sum(matin));
N1=1;
N2 = 20;
for i=N1:N2,
    for j=1:Nitr,
        [w,h]=NMF(matin',i);
        D = matin' - w*h;
        v(j,i-N1+1) = sum(sum(abs(D)));
    end
end
v = (min(v)/SS) * 100;
    