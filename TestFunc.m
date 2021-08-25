function TestFunc(func, flag)

	% load testing data
	fname = ['testData/' flag];
	load(fname);
	disp(fname);
	
	% test the functions
	if strcmp(flag,'a1') | strcmp(flag,'a2') 
		ResData = TestDataSet(func,'nLine',nLineSet,Data);
	end
	
	if strcmp(flag,'b1') | strcmp(flag,'b2')
		ResData = TestDataSet(func,'noise',noiseSet,Data);
	end
	
	if strcmp(flag,'c1') | strcmp(flag,'c2')
		ResData = TestDataSet(func,'outlier',vecRate,Data);
	end
	
	% store the experimental results
	fname = ['resData/' flag '_' func2str(func)];
	save(fname,'ResData');
    disp(fname);
    
function ResData = TestDataSet(func,str,vecID,Data)

	fprintf('\t\tmeanR\tmedianR\tmeanT\tmedianT\tcorrect-rate\ttime\n');
	
	nTest = length(Data{1});
	for i = 1:length(vecID)
		ID = vecID(i);
		fprintf('%s = %.2f   ', str, ID);
		
		% test
		t = tic;
		
		for j = 1:nTest
			D = Data{i}(j);
			[R_cw, T_cw] = func(D.p1, D.p2, D.P1_w, D.P2_w);
			ResData{i}.errR(j)= cal_rotation_err(R_cw, D.R_cw);
			ResData{i}.errT(j) = cal_translation_err(T_cw, D.T_cw);
		end
		
		t = toc(t);
	    ResData{i}.t = t / nTest * 1000;
	    
	    % show results
		printRes(ResData{i});
	end

function printRes(Res)

	fprintf('%f %f %f %f -- %.1f%% %.1fms\n', mean(Res.errR), ...
		median(Res.errR), ...
		mean(Res.errT), ...
		median(Res.errT), ...
		mean(Res.errR < 5 & Res.errT < 5) * 100, ...
		Res.t);
