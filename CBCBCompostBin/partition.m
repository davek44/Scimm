%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partition.m
%
% Partition a set of kmer count vectors into 2 sets based on
% the minimum normalized cut method used by CompostBin.
% Original CompostBin code is available at
% http://sites.google.com/site/souravc/compostbin
% and citable as "Sourav Chatterji, Ichitaro Yamazaki,
% Zhaojun Bai and Jonathan Eisen, CompostBin: A DNA
% composition-based algorithm for binning environmental
% shotgun reads , to appear in RECOMB 2008."
%
% Author: David Kelley
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

princ_comps = 3;

% load kmer counts
load kmers.dat;
% load indexes to consider
load indexes.dat;
% load constraints
%load constraints.dat;

% take only kmers of interest
mykmers = kmers(indexes,:);
%myconstraints = constraints(indexes,:);

num_neighbors = round(2 + .5*log(length(mykmers)));
max_nn = max(num_neighbors+3, round(1.25*num_neighbors));

% compute PCA
disp('Computing PCA');
[coeff,kmers_pc] = princomp(mykmers);

Wconn = false;
while ~Wconn && num_neighbors < max_nn
  % form nn graph and compute weights
  disp('Computing neighborhood');
  W = compW(kmers_pc(:,1:princ_comps), num_neighbors);
  % force symmetry
  W = max(W,W');

  % constrain reads
  %c1 = find(myconstraints == 0);
  %c2 = find(myconstraints == 1);
  %W(c1,c1) = 1.0;
  %W(c2,c2) = 1.0;
  %W(c1,c2) = 0.0;
  %W(c2,c1) = 0.0;

  if connected(W)
    Wconn = true;
  else
    num_neighbors = num_neighbors + 1;
    disp('Graph is disconnected, increasing number of neighbors');
  end
end

D = diag(sum(W));

% compute second smallest eigenval/vect
% 'sm' gave errors, esp. when the matrix has >1 conn comp
opt.disp = 0;
opt.tol = 1e-10;
opt.maxit = 1000;
disp('Computing Laplacian eigenvalues');
[V,E,FLAG] = eigs((D-W),D,2,'sa',opt);
if FLAG ~= 0
  disp('Eigenvector calculation did not converge');
  part = [0];
else
  % compute optimal partition
  disp('Computing optimal eigenvector partition');
  [part,Ncut] = split_optimal(W, V(:,2));
end

% print partition
fp = fopen('partition.txt','w');
for i=1:length(part)
  fprintf(fp, '%d\n', part(i));
end
fclose(fp);

% print cut
fp = fopen('ncut.txt','w');
fprintf(fp, '%f\n', Ncut(1,1));
fclose(fp);

exit;
