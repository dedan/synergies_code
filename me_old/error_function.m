
% I use this to compute the error. It is just the sum of the square 
% distance, from a prototype to all the vectors in the cluster

function error = error_function(proto, cluster)
error = 0;
for i = 1:size(cluster,1)
%   error = error + sum(proto - cluster(i,:))^2;
   error = error + pdist([proto; cluster(i,:)])^2;
end