% self-contained latex table generation for {table3} aka Table 2 in Sec 5.1,
% from saved MAT file data.  Writes to .tex inclusion file.
% Barnett 10/13/22.

%dir = [fileparts(mfilename('fullpath')) '/data'];   % local code writes to
dir = [fileparts(mfilename('fullpath')) '/efgp_table_results'];  % zip file

%out = 'table3_efgp_res_alex.tex';
out = 'table3_efgp_res.tex';

f=fopen(out,'w');

% table header
fprintf(f,'\\resizebox{\\textwidth}{!}{%%\n');
fprintf(f,"\\begin{tabular}{ccccccccccccc}\n");      % escaped "\" -> \\
fprintf(f,"   $N$ & $d$ & kernel & $\\varepsilon$ & $m$ & pre (s) & solve (s) & mean (s) & tot (s)  & iters & $\\textrm{EEPM}_{\\textrm{new}}$ & RMSE & RMSE$_{\\textrm{ex}}$\\\\ \n");
fprintf(f,"\\hline\n");

for ker = [1 2]
  for dim = 1:3
    if ker==1       
      fnam = 'se';  % for file
      knam = 'SE';  % for latex
      opts.tol = 1e-4;   % copied from se.m, should have been in info but wasnt
      if dim==3, opts.tol = 1e-3; end
    elseif ker==2
      fnam = 'mat';      % for file
      knam = 'Mat $1/2$';  % for latex
      if dim == 1        % copied from matern.m
        opts.tol = 1e-4;
      elseif dim == 2
        opts.tol = 1e-3;
      else
        opts.tol = 0.5*1e-2;  % philip's settled value
      end
    end
    for i=3:7
      N = 10^(i);
      % load info 
      filename = sprintf('%s_%gd_info_1e%g.mat', fnam, dim, i);
      load(fullfile(dir, filename));
      if ~isfield(info,'opts'), info.opts = opts; end   % for when appear
      if ~isfield(info,'opts_ref'), info.opts_ref = opts_ref; end
      % load errors
      filename = sprintf('%s_%gd_rms_err_1e%g.mat', fnam, dim, log10(N));
      load(fullfile(dir, filename));
      filename = sprintf('%s_%gd_eepm_err_1e%g.mat', fnam, dim, log10(N));
      load(fullfile(dir, filename));
      clear err_rmse_ex err_rms2    % deal with wrong var name
      filename = sprintf('%s_%gd_err_rmse_ex_1e%g.mat', fnam, dim, log10(N));
      load(fullfile(dir, filename));
      if ker==2 && ~exist('err_rmse_ex')
        err_rmse_ex = err_rms2;   % hack the wrong var name in last file
      end
      
      m = (numel(info.xis) -  1) / 2;
      % print for table

      % notice latex scientific notation formatting...
      fprintf(f,"$10^{%d}$ &$%d$ &%s &$%s$ &$%d$ & $%.3f$ & $ %.3f $ & $ %.3f $ & $ %.3f $ & $%d$ & $%s$ & $%s$ & $%s$ \\\\ \n", ...
              log10(N), dim, knam, scifmt(info.opts.tol,1), m, info.cpu_time.precomp, info.cpu_time.cg, info.cpu_time.mean, info.cpu_time.total, info.iter, scifmt(err_eepm), scifmt(err_rms), scifmt(err_rmse_ex));
    end
    fprintf(f,"\\hline \\hline \n");
  end
end
fprintf(f,"\\end{tabular}\n");
fprintf(f,"}%%");

% on alex's system, overwrites paper file...
system(sprintf('cp %s ../../../equispaced_fourier_gps/',out));
