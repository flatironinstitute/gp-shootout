A = load('data_efgp_get_xis_compare.mat');
B = load('data_efgp_get_xis_compare_cgtol_100.mat');



cc = ['b','r'];
mm = ['s','d'];
kernels = ["SE","Matern12","Matern32","Matern52"];
lstr = ["l=0.1","l=0.3"];
for dim=[1 2 3]
    if(dim == 1), tolvec = [1e-4; 1e-6]; end
    if(dim == 2), tolvec = [1e-3; 1e-5]; end
    if(dim == 3), tolvec = [1e-2; 1e-3]; end
    figure(dim);
    clf
    for ll=[1 2]
        for iker=[1 2 3 4]       
            isp = (ll-1)*4 + iker;
            subplot(2,4,isp)
            for itol=[1 2]
                
                avec = A.errs(:,iker,ll,itol,dim);
                bvec = B.errs(:,iker,ll,itol,dim);
                struse = [cc(itol) mm(1)];
                semilogy(avec,struse,'MarkerSize',5,'MarkerFaceColor',cc(itol)); hold on;
                struse = [cc(itol) mm(2)];
                semilogy(bvec,struse,'MarkerSize',5,'MarkerFaceColor',cc(itol));
                struse = [cc(itol) '--'];
                semilogy(1:3,ones(1,3)*tolvec(itol),struse);
   
            end
            struse = join([kernels(iker) ',  ' lstr(ll)]);
            title(struse);     
        end   
    end
    figure(dim);
    sgtitle(['dim=' int2str(dim)]);
end
%errs(imethod,iker,ll,itol,dim)