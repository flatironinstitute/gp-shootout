% self-contained latex table generation for {table3} aka Table 2 in Sec 5.1,
% from saved MAT file data.  Writes to .tex inclusion file.
% Barnett 10/13/22.

%dir = [fileparts(mfilename('fullpath')) '/data'];   % code writes to
dir = [fileparts(mfilename('fullpath')) '/results'];  % zip file

itab = 1



iNvals = [1 2 3 4 5];
if(itab ==1 )
    out = 'big_table_fix_sig.tex';
    isig = ones(size(iNvals)); 
else
    out = 'big_table_fix_noversig.tex';
    isig = iNvals;
end
f=fopen(out,'w');

% table header
fprintf(f,'\\resizebox{\\textwidth}{!}{%%\n');
fprintf(f,"\\begin{tabular}{ccccccccccc}\n");      % escaped "\" -> \\
fprintf(f,"  Alg & $\\sigma$ & $\\varepsilon$ & $N$ & $m$ & pre (s) & solve (s) & mean (s) & tot (s)  & iters & $\\textrm{EEPM}_{\\textrm{new}}$ & RMSE & RMSE$_{\\textrm{ex}}$\\\\ \n");
fprintf(f,"\\hline\n");

for i = 1:5
  fname3 = ['iNvals' int2str(iNvals(i)) '_isig' int2str(isig(i)) '_iNtrg3.mat'];
  iexp = -5
  epstr5 = sprintf('10^{%d}',iexp);
  iexp = -7
  epstr7 = sprintf('10^{%d}',iexp);
  if(i == 1 || i == 3)
    Nstr = sprintf('%.0f\\times 10^{%d}',A.N/10^(6+(i-1)/2), 6+(i-1)/2);   % escaped \
  else
    Nstr = sprintf('10^{%d}',log10(A.N));
  end
  fprintf(f,"$10^{%d}$ &$%d$ &%s &$%s$ &$%d$ & $%.3f$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $%d$ & $%s$ & $%s$ & $%s$ \\\\ \n", ...
              log10(N), dim, knam, epstr, m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, scifmt(err_eepm), scifmt(err_rms), scifmt(err_rmse_ex));

  fprintf(f,"\\hline \\hline \n");
end
fprintf(f,"\\end{tabular}\n");
fprintf(f,"}%%");

% overwrite paper file...
system(sprintf('cp %s ../../../equispaced_fourier_gps/',out));
