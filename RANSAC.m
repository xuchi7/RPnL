function [R, T] = RANSAC(p1, p2, P1_w, P2_w, nSAC)
% the line is expressed by the start and end points
% inputs:
%	 p1: 2d projection of the start point
%	 p2: 2d projection of the end point
%	 P1_w: 3d coordinates of the start point in world frame
%	 P2_w: 3d coordinates of the end point in world frame
% outputs:
%	 R: estimated rotation
%	 T: estimated translation

	if nSAC ~= 3 & nSAC ~= 4
		error 'n = 3 or 4';
	end

	n = length(p1);	
	
	nl = getProjNorm(p1,p2);

	% cal M
	M1 = kron([1 1 1], [nl nl]');
	M2 = kron([P1_w P2_w]', [1 1 1]);
	M = [M1 .* M2 [nl nl]'];
	
	thErr = 0.04;
	thN = max(0.4 * n, 3);
	maxN = -1;
	nIter = 0;
	while true
		nIter = nIter + 1;
		if nIter > 80
			if thN < 3
				break; 
			end
			thN = thN * 0.5;
			nIter = 1;
		end
				
		% get subset
		idx = getsubset(n,nSAC);
		sp1 = p1(:,idx);
		sp2 = p2(:,idx);
		sP1_w = P1_w(:,idx);
		sP2_w = P2_w(:,idx);
		
		if nSAC == 3
			sol = P3L(sp1,sp2,sP1_w,sP2_w);
		else % nSAC == 4
			[sol(1).R sol(1).T] = ASPnL(sp1,sp2,sP1_w,sP2_w);
		end
		
		% check and refine
		for j = 1:length(sol)						
			x = [sol(j).R(:); sol(j).T(:)];			
			Err = abs(M * x);
			idx = Err < thErr;
			idx = reshape(idx,[],2);
	    	idx = idx(:,1) & idx(:,2);
	    	resN = sum(idx);
	    	if resN > maxN
				R = sol(j).R;
				T = sol(j).T;
				maxN = resN;
			end
			if maxN > thN
				return;
			end
		end
	
	end

function subset = getsubset(n,m)

	subset = randi([1 n],m,1);
	subset = unique(subset);
	if length(subset) == m
		return;
	else
		subset = getsubset(n,m);
	end
