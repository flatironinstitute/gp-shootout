% self-contained latex table generation for {table3} aka Table 2 in Sec 5.1,
% from saved MAT file data.  Writes to .tex inclusion file.
% Barnett 10/13/22.
clear
%dir = [fileparts(mfilename('fullpath')) '/data'];   % code writes to
dir = [fileparts(mfilename('fullpath')) '/results'];  % zip file

itab = 2;



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
fprintf(f,"  Alg & $\\sigma$ & $\\varepsilon$ & $N$ & $m$ & iters & tot (s) & mem (GB) & $\\textrm{EEPM}$ & $\\textrm{EEPM}_{\\textrm{new}}$ & RMSE \\\\ \n");
fprintf(f,"\\hline\n");

for i = 1:5
  fname = ['iNvals' int2str(iNvals(i)) '_isig' int2str(isig(i)) '_iNtrg3.mat'];
  A = load(fname);
  iexp = -5;
  epstr5 = sprintf('10^{%d}',iexp);
  iexp = -7;
  epstr7 = sprintf('10^{%d}',iexp);
  if(i == 1 || i == 3)
    Nstr = sprintf('%.0f\\times 10^{%d}',A.N/10^(6+(i-1)/2), 6+(i-1)/2);   % escaped \
  else
    Nstr = sprintf('10^{%d}',log10(A.N));
  end
  
  
 fprintf(f,"%s & $%.2f$ & $%s$ & $%s$ & $%d$ & $%d$ & $%.0f$ & $%.1f$ & $%s$ & $%s$ & $%s$\\\\ \n", ...
              "EFGP", A.sig, epstr5, Nstr, (A.nxis_EFGP_5-1)/2, A.niter_EFGP_5, ...
          A.time_EFGP_5, A.mem_EFGP_5, scifmt(A.eepm_EFGP_5), ...
          scifmt(A.eepm_new_EFGP_5),scifmt(A.rms_EFGP_5));
 fprintf(f,"%s & $%.2f$ & $%s$ & $%s$ & $%d$ & $%d$ & $%.0f$ & $%.1f$ & $%s$ & $%s$ & $%s$\\\\ \n", ...
              "EFGP", A.sig, epstr7, Nstr, (A.nxis_EFGP_7-1)/2, A.niter_EFGP_7, ...
          A.time_EFGP_7, A.mem_EFGP_7, scifmt(A.eepm_EFGP_7), ...
          scifmt(A.eepm_new_EFGP_7),scifmt(A.rms_EFGP_7));
  if( i<= 3)
    fprintf(f,"%s & $%.2f$ & $%s$ & $%s$ & & & $%.0f$ & $%.1f$ & $%s$ & $%s$ & $%s$\\\\ \n", ...
              "FLAM", A.sig, epstr5, Nstr, ...
          A.time_FLAM_5, A.mem_FLAM_5, scifmt(A.eepm_FLAM_5), ...
          scifmt(A.eepm_new_FLAM_5),scifmt(A.rms_FLAM_5));
   fprintf(f,"%s & $%.2f$ & $%s$ & $%s$ & & & $%.0f$ & $%.1f$ & $%s$ & $%s$ & $%s$\\\\ \n", ...
              "FLAM", A.sig, epstr7, Nstr, ...
          A.time_FLAM_7, A.mem_FLAM_7, scifmt(A.eepm_FLAM_7), ...
          scifmt(A.eepm_new_FLAM_7),scifmt(A.rms_FLAM_7));
   
      
  end
  fprintf(f,"\\hline \\hline \n");
end
fprintf(f,"\\end{tabular}\n");
fprintf(f,"}%%");
% overwrite paper file...
system(sprintf('cp %s ../../../equispaced_fourier_gps/',out));
