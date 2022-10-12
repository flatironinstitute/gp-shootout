N = 3e3;        % problem size (small, matching naive, for now)
l = 0.1;        % SE kernel scale rel to domain [0,1]^dim, ie hardness of prob
sigma = 0.3;    % used to regress
sigmadata = sigma;   % meas noise, consistent case
freqdata = 3.0;   % how oscillatory underlying func? freq >> 0.3/l misspecified

opts = cell(3,1);
opts{1} = [];
opts{2} = [];
opts{3} = [];


opts{1}.use_integral = true;
opts{1}.cg_tol_fac = 100.0;

opts{2}.use_integral = false;
opts{2}.l2scaled = false;
opts{2}.cg_tol_fac = 100.0;


opts{3}.use_integral = false;
opts{3}.l2scaled = true;
opts{3}.cg_tol_fac = 100.0;
L = 1.0; shift = 200;   % arbitary, tests correct centering and L-box rescale
lvec = [0.1;0.3];


errs = zeros(3,4,2,2,3);
nxis = zeros(3,4,2,2,3);
iters = zeros(3,4,2,2,3);
ttime= zeros(3,4,2,2,3);
for dim = 1:3   % ..........
  unitvec = randn(dim,1); unitvec = unitvec/norm(unitvec);
  wavevec = freqdata*unitvec;    % col vec
  f = @(x) cos(2*pi*x'*wavevec + 1.3);   % underlying func, must give col vec
  rng(1); % set seed
  [x, meas, truemeas] = get_randdata(dim, N, f, sigmadata);    % x in [0,1]^dim
  x = L*x + (2*rand(dim,1)-1)*shift;                           % scale & shift
  if(dim == 1), tolvec = [1e-4; 1e-6]; end
  if(dim == 2), tolvec = [1e-3; 1e-5]; end
  if(dim == 3), tolvec = [1e-2; 1e-3]; end
     
  for itol = [1 2]
      fprintf('starting dim=%d,   itol=%d\n',dim,itol);
      tol = tolvec(itol);
      
  
      for ll = [1 2]
        l = lvec(ll);
        for iker = [1 2 3 4]
            if(iker == 1), ker = SE_ker(dim,L*l); end
            if(iker == 2), ker = Matern_ker(dim,0.5,L*l); end
            if(iker == 3), ker = Matern_ker(dim,1.5,L*l); end
            if(iker == 4), ker = Matern_ker(dim,2.5,L*l); end
            [ytrue, ytrg, ~] = naive_gp(x, meas, sigma^2, ker, [], opts{1});
            for imethod = [1 2 3]
                opts{imethod}.tol = tol;
                [y, ~, info] = EFGP(x, meas, sigma^2, ker, [], opts{imethod});
                iters(imethod,iker,ll,itol,dim) = info.iter;
                nxis(imethod,iker,ll,itol,dim) = (numel(info.xis)-1)/2;
                errs(imethod,iker,ll,itol,dim) = rms(y.mean-ytrue.mean);
                ttime(imethod,iker,ll,itol,dim) = info.cpu_time.total;
            end
        end
      end
  end
end

save('data_efgp_get_xis_compare_cgtol_100.mat','iters','nxis','errs','ttime');
