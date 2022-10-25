%% Define family of functions
clear
close('all')
x = chebfun('x',[0,1]);
logx = chebfun('log(x)',[0,1],'splitting','on');


lj = 1:0.04:40;
A = chebfun(exp(lj.*logx)); % x.^lj
[~,nfun] = size(A);

rng(1);
%% Create orthogonal basis for family of functions
B = randn(nfun);

% mix all functions, avoids the need to do pivoted QR
% Orthogonalize family of functions
[Q,R] = qr(A*B);

%% Use diagonal entries of R as an estimate for rank of family of functions
% Plot diagonal entries of R
figure(1)
clf
semilogy(abs(diag(R)),'k.');

%% Plot the first 6 orthogonal polynomials
figure(2)
clf
plot(Q(:,1:6),'LineWidth',1.6);

% Plot the same 6 orthogonal polynomials but with a log scale in x
% to see the structure close to the origin
figure(3)
clf
semilogx(Q(:,1:6),'LineWidth',1.6);


%% find number of basis functions which would approximate the full
% family of functions to given precision;

iind_min = min(find(abs(diag(R))<1e-9*abs(R(1,1)))); % don't try and use more,
                                           % than sqrt(eps_mach). Can run
                                           % into stability issues.
                                           % Best to do this in extended 
                                           % precision as a one time calc


% Find good interpolation nodes
iind = min(iind_min+1,nfun);

figure(4)
clf
plot(Q(:,iind+1),'LineWidth',1.6);
ff = Q(:,iind+1); hold on;
rr = roots(ff);
plot(rr,ff(rr),'k.');


ffuns = Q(:,1:iind);
B = ffuns(rr);

%% Test accuracy on random function in family 


lj_test = 1.1145;
fex = exp(lj_test*logx);

% Construct function samples
fvals = fex(rr);

% Apply values -> coefs matrix, note that inv(B) can be precomputed
% and stored since it is small
coefs = B\fvals;
finterp = ffuns*coefs;

figure(5)
plot(fex,'r-'); hold on; plot(finterp,'k-')

err_interp = norm(fex-finterp);
fprintf('error in interpolating an off-basis function=%d\n',err_interp);






