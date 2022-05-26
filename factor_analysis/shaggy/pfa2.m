function [] = pfa2(d,n)
% pfa.m   PRINCIPAL FACTOR ANALYSIS - a program designed to help determine
%         the number of significant factors in a data matrix.
%
% pfa(d) or pfa(d,n)
%
% d = data matrix
% n = number of principal factors to be saved in file temp.mat
%
% This program requires: slf.m and fabertbl.mat
% clc
% disp(' ')
% disp('                            PRINCIPAL FACTOR ANALYSIS')
% disp(' ')
% disp('                    @ copyright 2000 by Edmund R. Malinowski')
% disp(' ')
% disp('Type data identification statement:')
% disp(' ')
% id = input(' ','s');
% if isempty(id), id = (' '); end
id = (' ');
format short e
if nargin == 1, n = 0; end
[r,c] = size(d);
  sm = c; lg = r;
  if r < c, sm = r; lg = c; end
[u,s,v] = svd(d,0);
  for j = 1:sm
      ev(j) = s(j,j) * s(j,j);
      rev(j) = ev(j) / ((r-j+1)*(c-j+1));
  end
          sev(sm+1) = 0;
          sdf(sm+1) = 0;
      for k = sm:-1:2
          sev(k) = sev(k+1) + ev(k);
          sdf(k) = sdf(k+1) + (r-k+1) * (c-k+1);
      end
          for l = 1:sm-1
              m(l) = l;
              re(l) = sqrt(sev(l+1) / (lg * (sm-l)));
              ie(l) = re(l) * sqrt(l/sm);
              xe(l) = re(l) * sqrt((sm-l)/sm);
              ind(l) = re(l) / (sm-l)^2;
          end
m(sm) = sm; re(sm) = NaN; ind(sm) = NaN;
out(:,1) = m';
out(:,2) = ev';
out(:,3) = re';
out(:,4) = ind';
% Estimate the percent significance level using Malinowski's method
  for j = 1:sm-1
    f = (sdf(j+1) * ev(j)) / ((r-j+1) * (c-j+1) * sev(j+1));
    Msl(j) = slf(f,1,sm-j);
    if Msl(j) < 1e-2, Msl(j) = 0; end
  end
  Msl(sm) = NaN;
  out(:,5) = Msl';
% Estimate the percent significance level using Faber's method
for k = 1:sm; Fsl(k) = NaN; end
if lg <= 2000 & sm <= 100  
  for i=1:sm
    lgdf(i) = lg-i+1;
    smdf(i) = sm-i+1;
  end
  load fabertbl
  for k=1:sm-1

% df1(k) = table2(fabertbl,lgdf(k),smdf(k));

    df1(k) = interp2(fabertbl,lgdf(k),smdf(k));
    df2(k) = lgdf(k)*smdf(k) - df1(k);
  end
  for j=1:sm-1
    sev = sum(ev(j+1:sm));
    f(j) = (ev(j)/sev)*(df2(j)/df1(j));
    Fsl(j) =slf(f(j),df1(j),df2(j));
    if Fsl(j) < 1e-2, Fsl(j) = 0; end
   end
 end
  Fsl(sm) = NaN;
out(:,6) = Fsl';
% clc
disp('                          Principal Factor Analysis')
disp(id)
disp('       n           EV           RE           IND    %SL(Malinowski) %SL(Faber)')
disp(out)
% pause
% clc
save errors ev re ie xe ind Msl Fsl
% disp('The following results are stored in a file labelled "errors.mat":')
% disp(' ')
% disp(' ev = eigenvalues')
% disp(' re = real errors')
% disp(' ie = imbedded errors')
% disp(' xe = extracted errors')
% disp('ind = factor indicator function')
% disp('Msl = Malinowski`s percent significance level')
% disp('Fsl = Faber`s percent significance level')
% disp(' ')
% disp('To access the file type "load errors".')
if n > 0
  u = u(:,1:n);
  s = s(1:n,1:n);
  v = v(:,1:n);
  dr = u * s * v';
  de = dr - d;
  save temp u s v dr de
  disp(' ')
  disp('The following results are stored in a file labelled "temp.mat":')
  disp(' ')
  disp(' u = normalized abstract row factor matrix (each column is a factor)')
  disp(' v = normalized abstract column factor matrix (each column is a factor)')
  disp(' s = diagonal matrix of singular values (square roots of the eigenvalues)')
  disp(['dr = reproduced data matrix based on ',num2str(n),' principal factors'])
  disp('de = errors in data reproduction')
  disp(' ')
  disp('To access the file, type "load temp".')
end 

