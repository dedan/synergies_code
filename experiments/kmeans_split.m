function [centers, options, post,errlog] = ...
    kmeans_split(data, options, ncenters, metric_func)
%
% function [centers, options, post,errlog] = ...
%    kmeans_split(data, options, ncenters, metric_func)
%
%
% A hirarchic shell function for KMEANS function, 
%
% options
% 1: output
% 

[ndata, data_dim] = size(data);
if ~exist('metric_func')
  metric_func = 'dist2';
end


num_splits = log2(ncenters);
if options(1)>=0
    disp(num_splits);
end
if num_splits~=round(num_splits)
  ncenters = 2^round(num_splits)
  [centers, options, post,errlog] = ...
      kmeans_split(data, options, ncenters, metric_func);
  return;
end
% Initialize to one center
centers = mean(data);
post = ones(size(data,1),1);
delta = std(data)*.1;

for i=1:num_splits
  % Do the split
  if(options(1)>=0) disp(sprintf('split_iteration %d',i));end; 
  curr_s = size(centers,1);
  for j=1:size(centers,1)
    centers(j+curr_s,:) = centers(j,:) + delta.*rand(size(delta)); 
  end
  [centers,options,post,errlog] = kmeans(data,options,centers,metric_func);
end
  
