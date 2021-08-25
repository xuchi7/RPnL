function plotFigure(flag)

	S.titleSize = 16;
	S.winSize = [375 350];
	S.ColorList= {'r','b','m','k','c','g','y'};
	S.MarkerList= {'+','d','^','+','d','^','d'};
	S.styleidx = 1:7;
	vecFunc = {@Ansar,@Mirzaei,@ASPnL,@LPnL_DLT,@LPnL_Bar_LS,@LPnL_Bar_ENull};

	vecYLabel = {'Mean Rotation Error(deg)','Median Rotation Error(deg)',...
		'Mean Translation Error(%)','Median Translation Error(%)',...
		'Correct Rate(%)','Time(ms)'};
	
	vecTitle = cell(size(vecYLabel));
	timeRng = [0 300];
	rateRng = [0 100];

	if strcmp(flag,'a1') | strcmp(flag,'a2')
		X = 4:20;
		if strcmp(flag,'a1')
			rotaRng = [0 4];
			tranRng = [0 5];
			vecTitle = vecYLabel;
		elseif strcmp(flag,'a2')
			rotaRng = [0 40];
			tranRng = [0 50];
		end
		XLabel = 'Number of Lines';
	end
	
	if strcmp(flag,'b1') | strcmp(flag,'b2')
		X = 1:15;
		rotaRng = [0 20];
		tranRng = [0 50];
		XLabel = 'Noise(pixels)';
	end
	
	if strcmp(flag,'c1') | strcmp(flag,'c2')
		vecFunc = {@Mirzaei,@ASPnL,@RLPnL_LS,@RLPnL_ENull,@RANSAC3,@RANSAC4};
		S.MarkerList{1} = '^';
		S.MarkerList{4} = '.';
		S.styleidx = [2 3 5 6 1 4];
		X = 0.05:0.05:0.6;
		X = X * 100;
		rotaRng = [0 80];
		tranRng = [0 80];
		if strcmp(flag,'c1')
			rotaRng = [0 50];
			tranRng = [0 50];
		end
		XLabel = 'Outlier Rate(%)';
	end
	
	% load data
	for i = 1:length(vecFunc)
		 filename = [flag '_' func2str(vecFunc{i})];
		 Res(:,:,i) = loadData(['resData/' filename]);
		 S.MethodNames{i} = strrep(func2str(vecFunc{i}),'_','\_');
	end
	
	% draw figure
	vecYRng = {rotaRng, rotaRng, tranRng, tranRng, rateRng, timeRng};
	for i = 1:length(vecYRng)
		vecY = squeeze(Res(:,i,:));
		xPlot(i,X,vecY,vecTitle{i},XLabel,vecYLabel{i},S);
		ylim(vecYRng{i});
		if i == 5
			h=legend;
			set(h,'Loc','SouthEast');
		end
	end

function Res = loadData(filename)

	load(filename);
	n = length(ResData);
	for i = 1:n
		s = ResData{i};
		meanR(i) = mean(s.errR);
		medianR(i) = median(s.errR);
		meanT(i) = mean(s.errT) * 100;
		medianT(i) = median(s.errT) * 100;
		idx = s.errR < 5 & s.errT < 5;
		rate(i) = mean(idx) * 100;
		t(i) = s.t;
	end

	Res = [meanR(:) medianR(:) meanT(:) medianT(:) rate(:) t(:)];

function xPlot(fid,X,vecY,Title,XLable,YLable,S)

	% plot
	figure(fid); clf;

	plist= plot(X,vecY);
	for i= 1:length(plist)
		j = S.styleidx(i);
		set(plist(i),'Marker',S.MarkerList{j});
		set(plist(i),'Color',S.ColorList{j});
		set(plist(i),'MarkerFaceColor',S.ColorList{j});
	end
	legend(S.MethodNames);

	if ~isempty(Title)
		title(Title,'FontSize',S.titleSize);
	end

	xlabel(XLable);
	ylabel(YLable);
	xlim([X(1),X(end)]);

	set(gcf,'position',[100,100,S.winSize]);
	set(gcf,'PaperPositionMode','auto');

	return
