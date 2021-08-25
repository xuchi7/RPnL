function main()

% the experiments
	testCase a1; % centered (outlier-free)
	testCase a2; % uncentred (outlier-free)
	testCase b1; % nLine = 4
	testCase b2; % nLine = 5
	testCase c1; % centered (with outliers)
	testCase c2; % uncentred (with outliers)

function testCase(flag)
% flag: name of the experiment. We define 6 cases { a1, a2, b1, b2, c1, c2 }
% 	the meaning of the flags are as follows:
%		 a1; % centered (outlier-free)
%		 a2; % uncentred (outlier-free)
%		 b1; % nLine = 4
%		 b2; % nLine = 5
%		 c1; % centered (with outliers)
%		 c2; % uncentred (with outliers)
	
	% functions compared
	if strcmp(flag,'a1') | strcmp(flag,'a2') | strcmp(flag,'b1') | strcmp(flag,'b2') 
		funcnames = {@LPnL_DLT,@LPnL_Bar_LS,@LPnL_Bar_ENull,@Ansar,@Mirzaei,@ASPnL};
	elseif strcmp(flag,'c1') | strcmp(flag,'c2')
		funcnames = {@Mirzaei,@ASPnL,@RLPnL_LS,@RLPnL_ENull,@RANSAC3,@RANSAC4};
	else
		error 'the flag of the case is not defined!'
	end
	
	% generate testing data, and store it in "testData" folder
	genTestData(flag);
	% load testing data, test the functions, and store the results in "resData" folder
	disp(flag);
	for j = 1:length(funcnames)
		fprintf('\t%s\n',func2str(funcnames{j}));
		TestFunc(funcnames{j}, flag);
	end
	% read results in "resData" folder, and plot the figures
	plotFigure(flag);

