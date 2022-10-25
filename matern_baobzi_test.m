addpath /Users/mrachh/local/baobzi/share/baobzi/matlab/
dim = 1;
order = 6;

rho = 10;
nu = 3.8;

dmax = rho/sqrt(2*nu)*log(1.0/eps);

center = dmax/2;
hl = dmax;
tol = 1E-10;
minimum_leaf_fraction = 0.0;
split_multi_eval = 1;


func_approx = baobzi('new', 'matern_fun', dim, order, center, hl, tol, minimum_leaf_fraction, split_multi_eval);
display(func_approx.eval([0.1]))
func_approx.save('matern_rho10_nu3.8.baobzi');
clear func_approx;
func_approx = baobzi('restore', 'matern_rho10_nu3.8.baobzi');
func_approx.stats()
display(func_approx.eval([0.25, 0.25]))

n = 5e5;
x = rand([n, 1]);
tic; y1 = func_approx.eval(x); t1= toc;

prefac = 2^(1-nu)/gamma(nu);

% Setting nuber of threads to 1 since baobzi is single threaded
nth = maxNumCompThreads(1);
tic; 
xsc = sqrt(2*nu)/rho*x;
y2 = prefac*xsc.^(nu).*besselk(nu,xsc);
t2 = toc;
fprintf('-----------------------------------------------------\n');
fprintf('speed using baobzi Meval/core/sec=%0.1f\n',n/t1/1e6)
fprintf('time to evaluate matern kernel natively Meval/core/sec=%0.1f\n\n',n/t2/1e6)
error_mat = norm(y2-y1)/norm(y2);
fprintf('speed up factor=%0.1f\n\n',t2/t1);
fprintf('error in evaluation=%d\n',error_mat);

