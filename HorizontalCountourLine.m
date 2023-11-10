tic;


if ~exist(outputFolder_simulation,'dir')
    mkdir(outputFolder_simulation);
end

CurrentFig = casefoldername;


DPOSIT_List = nan(10000,1);
DPOSIT_temp = 0;

dispfd_List = nan(10000,1);
dispfd_temp = 0;

sizef_List = nan(10000,1);
sizef_temp = 0;

TIME_list = nan(10000,1);
TIME_temp = 0;

for FrameNumber = 100:100 %NumberOfFrames   % Finish loop once for each movie frame.

  
  cd(outputFolder_simulation);
  load('TopoINfO.mat');
  
  erosion = -0.01; % erosion between slices
  rotationDegree = 0; % rotation between slices
  
  maxHeight = max(max(max(ZBEDDraw)));
  minHeight = min(min(min(ZBEDDraw)));
  maxHeightStart = max(max(ZBEDDraw(:,:,1)));
  minHeightStart = min(min(ZBEDDraw(:,:,1)));
  maxHeightEnd = max(max(ZBEDDraw(:,:,end)));
  minHeightEnd = min(min(ZBEDDraw(:,:,end)));
  
  sliceN = floor(((minHeightEnd - minHeightStart)/(maxHeightStart-minHeightStart)-1.001)/erosion) - 1;

  % dTrend = rotationDegree : rotationDegree : rotationDegree * sliceN;
  dTrend = zeros(1, sliceN); % no rotation

  ContourLinesList = cell(1,sliceN);
  zrefList = NaN(1,sliceN);
  
  
  fg = figure(1);
  for slicei = 1 : sliceN
      
      DuneInit;
      
      % timerecorded = NaN(1, NBEDSH);
      contourTable = cell(1,NBEDSH);
      
      ZHORIZ = ZHORIZ + erosion*slicei;
      TRENDF = TRENDF + rotationDegree*slicei;
      TRENDS = TRENDS + rotationDegree*slicei;
      TRENDT = TRENDT + rotationDegree*slicei;
%       ZHORIZ = ZHORIZ + dZHO(slicei);
%       TRENDF = TRENDF + dTrend(slicei);
%       TRENDS = TRENDS + dTrend(slicei);
%       TRENDT = TRENDT + dTrend(slicei);

      
      for N = 1 : NBEDSH
        TIME = 1 + dT(FrameNumber) - N;
        TIMEEnd = dT(FrameNumber);
        TIMEStart = 1 + dT(FrameNumber) - NBEDSH;
        
        
        if MORVRT == false; break; end
        if FirstRun
          DuneTopo;
          % sandsurface;   

          if max(max(z)) > zref
             MORHRZ = true;
          end

          zcont = z;
          FirstRun = false;

        end


        if MORHRZ % horizontal slice (zref) cutting through modeled topography 
          x = XSurf;
          y = YSurf;
          DuneTopo;
          zcont = min(z,zcont); % could be eroded by later time

          maxz = max(max(zcont));
          minz = min(min(zcont));
          if mod(TIME,INTXBD) == 0 % INTXBD: time interval between the plotting of foresets
            
              set(0,'CurrentFigure',fg);
              [c,h]=contour3(x,y,zcont,[zref zref]);                  %surface contours
              
              hold on;
              xlim([min(min(XSurf)) max(max(XSurf))]);
              ylim([min(min(YSurf)) max(max(YSurf))]);
              zlim([minHeight maxHeight]);
              
              
            if ~isempty(c)
                contourTableThisSurface = getContourLineCoordinates(c); % level, group (line number), x, y
                contourTable{1,N} = contourTableThisSurface;
                
                % set(h,'linewidth', 1.25,'edgecolor',[0.81 0.63 0.42])    %surface contour colors  - change linewidth here!! 1.25 is thick line
            end

          end

          if max(max(zcont)) <= zref; MORHRZ = false; end
        end

      end  % End of loop to finish each image
      drawnow;
      
      EmptyCells = find(cellfun(@isempty,contourTable));
      contourTable(EmptyCells) = []; % delete empty cells

      tempGroupn = 0;
      linenumber = 1;
      AllContourLines = NaN(100000,4);
      for i = 1 : size(contourTable,2)
          thiscell = contourTable{i};
          groupnumber = max(thiscell.Group);
          cellmatrix = table2array(thiscell);
          cellmatrix(:,2) = cellmatrix(:,2) + tempGroupn;

          AllContourLines(linenumber: linenumber + size(cellmatrix,1)-1,:) = cellmatrix;

          tempGroupn = tempGroupn + groupnumber;
          linenumber = linenumber + size(cellmatrix,1);

      end
      [rdel,~] = find(isnan(AllContourLines));
      AllContourLines(rdel,:) = [];

      ContourLinesList{1,slicei} = AllContourLines;
      zrefList(1,slicei) = zref;


  
  hold off;
  end
  
  % SAVE THE MATRIX FOR HORIZONTAL SLICES
  cd(outputFolder_simulation);
  save('ContourLineINfO','zrefList','ContourLinesList');
  
  
  
 

%  % draw horizontal section
%  I = 1;
%  thisHorContLines = ContourLinesList{I};
%  for checki = 1 : max(thisHorContLines(:,2))
%      [rcont,~] = find(thisHorContLines(:,2)==checki);
%      thisline = thisHorContLines(rcont,:);
%      plot3(thisline(:,3),thisline(:,4),thisline(:,1));
%      hold on;
%  end

end  % End of loop to make movie sequence.
toc;
