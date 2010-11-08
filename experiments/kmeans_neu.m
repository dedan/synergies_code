function [centres, options, post, errlog] = kmeans(...
    data, options, centres, metric_func)
%
% KMEANS Trains a K-Means cluster model.
%
%	Description
%	 CENTRES = KMEANS(DATA, OPTIONS,CENTRES) uses the batch K-means
%	algorithm to set the centres of a cluster model. The matrix DATA
%	represents the data which is being clustered, with each row
%	corresponding to a vector. The sum of squares error function is used.
%	The point at which a local minimum is achieved is returned as
%	CENTRES.  The error value at that point is returned in OPTIONS(8).
%
%	[CENTRES, OPTIONS, POST, ERRLOG] = KMEANS(CENTRES, DATA, OPTIONS)
%	also returns the cluster number (in a one-of-N encoding) for each
%	data point in POST and a log of the error values after each cycle in
%	ERRLOG.    The optional parameters have the following
%	interpretations.
%
%	OPTIONS(1) is set to 1 to display error values; also logs error
%	values in the return argument ERRLOG. If OPTIONS(1) is set to 0, then
%	only warning messages are displayed.  If OPTIONS(1) is -1, then
%	nothing is displayed.
%
%	OPTIONS(2) is a measure of the absolute precision required for the
%	value of CENTRES at the solution.  If the absolute difference between
%	the values of CENTRES between two successive steps is less than
%	OPTIONS(2), then this condition is satisfied.
%
%	OPTIONS(3) is a measure of the precision required of the error
%	function at the solution.  If the absolute difference between the
%	error functions between two successive steps is less than OPTIONS(3),
%	then this condition is satisfied. Both this and the previous
%	condition must be satisfied for termination.
%
%	OPTIONS(14) is the maximum number of iterations; default 100.
%
%	See also
%	GMMINIT, GMMEM
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

init_centres = 0;
toplot = 0;

[ndata, data_dim] = size(data);
if (length(centres)~=1)
  [ncentres, dim] = size(centres);
else
  dim = data_dim;
  ncentres = centres;
  init_centres = 1;
end


if ~exist('metric_func')
  metric_func = 'dist2';
end

if (dim ~= data_dim)
  dim
  data_dim
  error('Data dimension does not match dimension of centres')
end

if (ncentres > ndata)
  error('More centres than data')
end

% Sort out the options
if (options(14))
  niters = options(14);
else
  niters = 100;
end

store = 0;
if (nargout > 3)
  store = 1;
  errlog = zeros(1, niters);
end

% Check if centres and posteriors need to be initialised from data
% Amir: Changed this a bit so centres are initialized to means of
% groups of data

if (options(5) == 1 | init_centres)
  % Do the initialisation
  clear centres;
%  perm = randperm(ndata);
%  centres = data(perm(1:ncentres),:);
  for i=1:ncentres
    perm = randperm(ndata);
    perm = perm(1:floor(ndata/ncentres));
    centres(i,:) = mean(data(perm,:));
  end
end

% Matrix to make unit vectors easy to construct
id = eye(ncentres);

% Main loop of algorithm
for n = 1:niters
  % Save old centres to check for termination
  old_centres = centres;
  
  % Calculate posteriors based on existing centres
  
  
%  d2 = dist2(data, centres);
  d2 = feval(metric_func,data,centres);
  
  % Assign each point to nearest centre
  [minvals, index] = min(d2', [], 1);
  post = id(index,:);

  num_points = sum(post, 1);
  % Adjust the centres based on new posteriors
  for j = 1:ncentres
    if (num_points(j) > 0)
      centres(j,:) = sum(data(find(post(:,j)),:), 1)/num_points(j);
    end
  end
  % Error value is total squared distance from cluster centres
  e = sum(minvals);
  if store
    errlog(n) = e;
  end
  if options(1) > 0
    fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
  end

  if n > 1
    % Test for termination
%    if max(max(abs(centres - old_centres))) < options(2) & ...
%        abs(old_e - e) < options(3)
     if abs(old_e - e) < options(3)
      options(8) = e;
      return;
    end
  end
  old_e = e;
  
  % Two dimensional plot of centres
  if (toplot)
    plot(centres(:,1),centres(:,2),'*');
    pause(.1);
  end
end

% If we get here, then we haven't terminated in the given number of 
% iterations.

options(8) = e;
%if (options(1) >= 0)
%  disp('Warning: Maximum number of iterations has been exceeded');
%end
