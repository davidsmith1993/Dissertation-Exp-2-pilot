function Y2 = genMsamples(u,K,n,catnum,m)
%function Y2 = genMsamples(u,K,n,catnum,m)
%
%  generates m categories of 'n' normally distributed data points.
%
%  Parameters:
%    u is a 2xm matrix containing the population mean
%    K is a 2x2 matrix containing the population covariance 
%      matrix (assumes equal K across categories)
%    n is the number of data points
%    catnum is a vector of category labels
%	  m is the number of categories (optional; defaults to gensample.m when omitted/m=1)
%
%  outputs a structure of size m; Row format of output:  [catnum x y ...]

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 16-April-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%	5/21/01			include algorithm to assure that sample stats match pop params

if nargin<5
   m=1;
end

dim = size(u,1);
Y2 = zeros(1,dim+1);

for i=1:m
	% Step 1: Generate X, a matrix containing N normally
	% distributed mean 0 covariance I random numbers.
	X = randn(dim,n(i));

	% Step 2: Transform X to have covariance K and mean u. 
	%This is accomplished via the linear transformation Y = AX+B
	Ksamp = cov(X');

	A=(K^.5)*(Ksamp^-.5);
	b=u(:,i)-A*mean(X')';
	Y=A*X+b*ones(1,n(i));


	% Check the covariance and mean of Y to see that it is close to
	% what we expected.
	Sample_K = cov(Y');
%    if Sample_K ~= K; 
%        pop_K=K
%        Sample_K = Sample_K
%        error('population and sample covariance matrices not equal'); 
%    end
	meanY = mean(Y');
%    if meanY' ~= u(:,i); 
%        pop_mean = u(:,i)
%        sample_mean = meanY
%        error('population and sample means not equal'); 
%    end
	
    % Also, add a row to Y to indicate the category from which the points
	% were generated.
   Y = [catnum(i)*ones(n(i),1) Y']; 
   Y2 = [Y2;Y];
end

%remove zeros from Y2
Y2=Y2(2:size(Y2,1),:);