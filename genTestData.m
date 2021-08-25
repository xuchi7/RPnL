function genTestData(flag)
	
	disp(['generating testing data - ' flag]);
	
	isUC = 0; % is uncentred or not
	nTest = 500; % number of tests
	noiseStd = 5; % noise level
	rateOutlier = 0; % out lier rate
	
	if strcmp(flag,'a1') | strcmp(flag,'a2')
		nLineSet = 4:20;
		if strcmp(flag,'a2')
			isUC = 1; % uncentred data
		end
		for i = 1:length(nLineSet)
			Data{i} = prepareTestData(nTest, nLineSet(i), noiseStd, rateOutlier, isUC);
		end
		% save testing data to "testData" folder
		fname = ['testData/' flag];
		save(fname, 'nTest', 'nLineSet', 'noiseStd', 'Data');
		whos('-file', fname);
	end

	if strcmp(flag,'b1') | strcmp(flag,'b2')
		nLine = 4;
		if strcmp(flag,'b2')
			nLine = 5;
		end
		noiseSet = 1:15;
		for i = 1:length(noiseSet)
			Data{i} = prepareTestData(nTest, nLine, noiseSet(i), rateOutlier, isUC);
		end
		% save testing data to "testData" folder
		fname = ['testData/' flag];
		save(fname, 'nTest', 'nLine', 'noiseSet', 'Data');
		whos('-file', fname);
	end
	
	if strcmp(flag,'c1') | strcmp(flag,'c2')
		if strcmp(flag,'c2')
			isUC = 1; % uncentred
		end
		nLine = 100;
		vecRate = 0.05:0.05:0.6; % outlier rate
		for i = 1:length(vecRate)
			Data{i} = prepareTestData(nTest, nLine, noiseStd, vecRate(i), isUC);
		end
		% save testing data to "testData" folder
		fname = ['testData/' flag];
		save(fname, 'nTest', 'nLine', 'noiseStd', 'Data', 'vecRate');
		whos('-file', fname);
	end		
	
function [pt,PT] = randPts(n,w,h,f,noiseStd,isUC)

	if isUC % uncentred
		u = rand(1,n) * w * 0.25;
		v = rand(1,n) * h * 0.25;
	else % centred
		u = rand(1,n) * w;
		v = rand(1,n) * h;
	end

	% depth range [4,8]
	depth = 4*rand(1,n)+4;
	
	% normalize
	u = (u - (w*0.5)) / f;
	v = (v - (h*0.5)) / f;

	% 2d
	pt = [u; v];

	% 3d
	PT = [u; v; ones(1,n)] .* kron(ones(3,1),depth);

	%add some noise to the endpoints
	if noiseStd > 0
		pt = pt + noiseStd*(rand(2,n)-0.5);
	end
		
function [i1,i2] = randOutlier(n,rateOutlier)

	m = round(n * rateOutlier);
	i0 = randperm(n);
	i1 = i0(1:m);
	i2 = i1(randperm(m));

function D = prepareTestData(nTest, nLine, noiseStd, rateOutlier, isUC)

	% parameters
    width = 640;
    height = 480;
	focal = 800;
	noiseStd   = noiseStd / focal; 

	% generate Data
	for kk= 1:nTest		    
		
		% lines
		[p1 P1_c] = randPts(nLine, width, height, focal, noiseStd, isUC);
		[p2 P2_c] = randPts(nLine, width, height, focal, noiseStd, isUC);
				
		%Rotation matrix from camera to world frame R_cw
		R_cw = randR();
		R_wc = R_cw.';
		
		% translation from camera to world frame T_cw
		T_cw = mean([P1_c P2_c],2);
		T_wc = -R_wc * T_cw;
				
		% world frame
		P1_w = R_wc * P1_c + kron(ones(1,nLine), T_wc);
		P2_w = R_wc * P2_c + kron(ones(1,nLine), T_wc);	

		% add outliers
		if rateOutlier > 0
			[i1,i2] = randOutlier(nLine,rateOutlier);
			p1(:,i1) = p1(:,i2);
			p2(:,i1) = p2(:,i2);
		end
		
		% save data
		D(kk).p1 = p1;
		D(kk).p2 = p2;
		D(kk).P1_w = P1_w;
		D(kk).P2_w = P2_w;
		D(kk).R_cw = R_cw;
		D(kk).T_cw = T_cw;
		
    end
        
function R = randR

    theta_x = 2*rand(1)-1;
    theta_y = 2*rand(1)-1;
    theta_z = 2*rand(1)-1;
    r_x = [1,0,0; 0, cos(theta_x), - sin(theta_x); 0, sin(theta_x), cos(theta_x) ];
    r_y = [cos(theta_y), 0, sin(theta_y); 0, 1, 0; - sin(theta_y), 0, cos(theta_y)];
    r_z = [cos(theta_z), - sin(theta_z), 0; sin(theta_z), cos(theta_z), 0; 0, 0, 1];
    R = r_z * r_y * r_x;
