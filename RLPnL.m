function [R, T] = RLPnL(p1, p2, P1_w, P2_w, isENull)
% the line is expressed by the start and end points
% inputs:
%	 p1: 2d projection of the start point
%	 p2: 2d projection of the end point
%	 P1_w: 3d coordinates of the start point in world frame
%	 P2_w: 3d coordinates of the end point in world frame
% outputs:
%	 R: estimated rotation
%	 T: estimated translation

	n = length(p1);
	if n < 8
		[R, T] = LPnL(p1, p2, P1_w, P2_w, ones(n,1) > 0, isENull);
		return;
	end

	% get base points
	ct_w = mean([P1_w P2_w],2);
	BP_w = [eye(3) zeros(3,1)] + kron([1 1 1 1],ct_w);
	alp1_w = getAlpha(BP_w,P1_w);
	alp2_w = getAlpha(BP_w,P2_w);

	% normal of porjected lines
	nl = getProjNorm(p1,p2);
  
    % matrix
    M1 = kron([1 1 1 1], [nl nl]');
    M2 = kron([alp1_w alp2_w]', [1 1 1]);
    M = M1 .* M2;
    D = ones(2*n,1);

	thErr = 0.04;
	Err0 = 100;
	for j = 1:5
		if isENull
			% ENull
			[R,T,Err1,vecErr] = solveENull(BP_w, M, D);
		else
			% LS
			[R,T,Err1,vecErr] = solveLS(BP_w, M, D);
		end
			
		% terminate
		if (Err1 > Err0)
			return;
		end
		% update Err0
		Err0 = Err1;
		
		% update D
		idx = vecErr < max(thErr,2*Err1);
	    idx = reshape(idx,[],2);
	    idx = idx(:,1) & idx(:,2);
	    D = [idx; idx];
	    m = sum(idx);
		if m > 0.9*n
			return;
		end		
		if m < max(6,n*0.2)
			return;
		end	
    end    
       
function val = getSortValue(vecerr, percent)
	vecerr_ = sort(vecerr);
	val = vecerr_(ceil(length(vecerr_)*percent));
    
function [R,T,Err1,vecErr] = solveLS(BP_w, M, D)

	MTM = M' * diag(D) * M;
	[XV XD] = xeig(MTM);

	X1 = XV(:,1);
	[X lm] = normBP(X1);
	vecErr = abs(M*X);
	Err1 = getSortValue(vecErr,0.3);
	
	BP_c = reshape(X,3,4);
	[R,T] = getRT(BP_w.',BP_c.');	
	
function [R,T,Err1,vecErr] = solveENull(BP_w, M, D)
% ====================== Similar to EPnP+GN solver ========================
global flagpath
if isempty(flagpath)
    addpath epnp;
    flagpath = 1;
end

%Compute kernel M
MTM = M' * diag(D) * M;
[XV XD] = xeig(MTM);
Km = XV(:,4:-1:1);

%1.-Solve assuming dim(ker(M))=1. X=[Km_end];------------------------------
X1=Km(:,end);
[X lm] = normBP(X1);
BP_c = reshape(X,3,4);
[R T] = getRT(BP_w.',BP_c.');
err(1)=x_calErr(BP_w,M,R,T,BP_c);

sol(1).R=R;
sol(1).T=T;
sol(1).error=err(1);
sol(1).betas=[1];
sol(1).sc=1/lm;
sol(1).Kernel=X1;
sol(1).BP_c = BP_c;

%2.-Solve assuming dim(ker(M))=2------------------------------------------
Km1=Km(:,end-1);
Km2=Km(:,end);

%control points distance constraint
D=compute_constraint_distance_2param_6eq_3unk(Km1,Km2);
dsq=define_distances_btw_control_points();
betas_=inv(D'*D)*D'*dsq;
beta1=sqrt(abs(betas_(1)));
beta2=sqrt(abs(betas_(3)))*sign(betas_(2))*sign(betas_(1));
X2=beta1*Km1+beta2*Km2;

[X lm] = normBP(X2);
BP_c = reshape(X,3,4);
[R T] = getRT(BP_w.',BP_c.');
err(2)=x_calErr(BP_w,M,R,T,BP_c);

sol(2).R=R;
sol(2).T=T;
sol(2).error=err(2);
sol(2).betas=[beta1,beta2];
sol(2).sc=1/lm;
sol(2).Kernel=[Km1,Km2];
sol(2).BP_c = BP_c;

%3.-Solve assuming dim(ker(M))=3------------------------------------------
Km1=Km(:,end-2);
Km2=Km(:,end-1);
Km3=Km(:,end);

%control points distance constraint
D=compute_constraint_distance_3param_6eq_6unk(Km1,Km2,Km3);
dsq=define_distances_btw_control_points();
betas_=inv(D)*dsq;
beta1=sqrt(abs(betas_(1)));
beta2=sqrt(abs(betas_(4)))*sign(betas_(2))*sign(betas_(1));
beta3=sqrt(abs(betas_(6)))*sign(betas_(3))*sign(betas_(1));

X3=beta1*Km1+beta2*Km2+beta3*Km3;

[X lm] = normBP(X3);
BP_c = reshape(X,3,4);
[R T] = getRT(BP_w.',BP_c.');
err(3)=x_calErr(BP_w,M,R,T,BP_c);

sol(3).R=R;
sol(3).T=T;
sol(3).error=err(3);
sol(3).betas=[beta1,beta2,beta3];
sol(3).sc=1/lm;
sol(3).Kernel=[Km1,Km2,Km3];
sol(3).BP_c = BP_c;

%4.-Solve assuming dim(ker(M))=4------------------------------------------
Km1=Km(:,end-3);
Km2=Km(:,end-2);
Km3=Km(:,end-1);
Km4=Km(:,end);

D=compute_constraint_distance_orthog_4param_9eq_10unk(Km1,Km2,Km3,Km4);
dsq=define_distances_btw_control_points();
lastcolumn=[-dsq',0,0,0]';
D_=[D,lastcolumn];
Kd=null(D_);

P=compute_permutation_constraint4(Kd);
lambdas_=kernel_noise(P,1);
lambda(1)=sqrt(abs(lambdas_(1)));
lambda(2)=sqrt(abs(lambdas_(6)))*sign(lambdas_(2))*sign(lambdas_(1));
lambda(3)=sqrt(abs(lambdas_(10)))*sign(lambdas_(3))*sign(lambdas_(1));
lambda(4)=sqrt(abs(lambdas_(13)))*sign(lambdas_(4))*sign(lambdas_(1));
lambda(5)=sqrt(abs(lambdas_(15)))*sign(lambdas_(5))*sign(lambdas_(1));

betass_=lambda(1)*Kd(:,1)+lambda(2)*Kd(:,2)+lambda(3)*Kd(:,3)+lambda(4)*Kd(:,4)+lambda(5)*Kd(:,5);
beta1=sqrt(abs(betass_(1)));
beta2=sqrt(abs(betass_(5)))*sign(betass_(2));
beta3=sqrt(abs(betass_(8)))*sign(betass_(3));
beta4=sqrt(abs(betass_(10)))*sign(betass_(4));
X4=beta1*Km1+beta2*Km2+beta3*Km3+beta4*Km4;

[X lm] = normBP(X4);
BP_c = reshape(X,3,4);
[R T] = getRT(BP_w.',BP_c.');
err(4)=x_calErr(BP_w,M,R,T,BP_c);

sol(4).R=R;
sol(4).T=T;
sol(4).error=err(4);
sol(4).betas=[beta1,beta2,beta3,beta4];
sol(4).sc=1/lm;
sol(4).Kernel=[Km1,Km2,Km3,Km4];
sol(4).BP_c = BP_c;
    
%5.-Gauss Newton Optimization------------------------------------------------------ 
[min_err,best_solution]=min(err);
R=sol(best_solution).R;
T=sol(best_solution).T;
Betas=sol(best_solution).betas;
sc=sol(best_solution).sc;
 
if best_solution==1
    Betas=[0,0,0,Betas];
elseif best_solution==2
    Betas=[0,0,Betas];
elseif best_solution==3
    Betas=[0,Betas];
end

Km1=Km(:,end-3);
Km2=Km(:,end-2);
Km3=Km(:,end-1);
Km4=Km(:,end);
Kernel=[Km1,Km2,Km3,Km4];

%refine the solution iterating over the betas
Beta0=Betas/sc;
[R_opt,T_opt,BP_c_opt,err_opt,iter]=x_gauss_newton(Kernel,BP_w,Beta0,M);

%Just update R,T,Xc if Gauss Newton improves results (which is almost
%always)
if err_opt<min_err    
    R=R_opt;
    T=T_opt;
    BP_c = BP_c_opt;
end

[Err1 vecErr]=x_calErr(BP_w,M,R,T,BP_c);

opt.Beta0=Beta0;
opt.Kernel=Kernel;
opt.iter=iter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Err vecErr] = x_calErr(BP_w,M,R,T,BP_c)

if (1)
    BP_c = R * BP_w + kron([1 1 1 1],T);    
end
	vecErr = abs(M*BP_c(:));
    Err = getSortValue(vecErr,0.3);
	
function [R_opt,T_opt,BP_c_opt,err_opt,iter]=x_gauss_newton(Km,BP_w,Beta0,M)

    n=size(Beta0,2);

    [Beta_opt,err,iter]=gauss_newton(Km,BP_w.',Beta0);

    X=zeros(12,1);
    for i=1:n
       X=X+Beta_opt(i)*Km(:,i); 
    end

    [X lm] = normBP(X);
    BP_c_opt = reshape(X,3,4);
    [R_opt,T_opt] = getRT(BP_w.',BP_c_opt.');
    err_opt=x_calErr(BP_w,M,R_opt,T_opt,BP_c_opt);
