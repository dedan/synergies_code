function [W,H,score] = my_svd(V,k)

[W,D,H] =  svds(V,k);
H= D *H';

dif = V - W*H; 
 score = sum(sum(sqrt(dif.*dif)));

end


