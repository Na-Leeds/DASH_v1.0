% get data for horizontal slice
if ~exist(outputFolder_simulation,'dir')
    mkdir(outputFolder_simulation);
end
NumberOfFrames = 300;
GridSize = 100;
SurfGridSpace = 0.5;
EdgeGridSpace = 0.5;

CurrentFig = casefoldername;

DPOSIT_List = nan(10000,1);
DPOSIT_temp = 0;

dispfd_List = nan(10000,1);
dispfd_temp = 0;

sizef_List = nan(10000,1);
sizef_temp = 0;

TIME_list = nan(10000,1);
TIME_temp = 0;


for FrameNumber = NumberOfFrames: NumberOfFrames

  
%   erosion = -0.025; % erosion between slices
%   rotationDegree = 0; % rotation between slices
%   sliceN = 100; % number of horizontal slices to draw
%   
%   dZHO = erosion : erosion : erosion * sliceN;
%   % dTrend = rotationDegree : rotationDegree : rotationDegree * sliceN;
%   dTrend = zeros(1, sliceN); % no rotation
% 
%   ContourLinesList = cell(1,sliceN);
%   zrefList = NaN(1,sliceN);
  
  
  fg = figure(1);

      
  DuneInit;

  
  ZBEDDraw = NaN(gridnumber,gridnumber,NBEDSH);
  apparentDipDraw = NaN(gridnumber,gridnumber,NBEDSH);
  apparentDipSignDraw = NaN(gridnumber,gridnumber,NBEDSH);
  ZErosion = NaN(gridnumber,gridnumber,NBEDSH); 
  timerecorded = NaN(1, NBEDSH);

%   ZHORIZ = ZHORIZ + erosion*slicei;
%   TRENDF = TRENDF + rotationDegree*slicei;
%   TRENDS = TRENDS + rotationDegree*slicei;
%   TRENDT = TRENDT + rotationDegree*slicei;
  
%       ZHORIZ = ZHORIZ + dZHO(slicei);
%       TRENDF = TRENDF + dTrend(slicei);
%       TRENDS = TRENDS + dTrend(slicei);
%       TRENDT = TRENDT + dTrend(slicei);


  temp = 0;
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


    if MORVRT
          x = XSurf;
          y = YSurf;
          DuneTopo;
          zcont = min(z,zcont); % could be eroded by later time
          
          
          
          
          if mod(TIME,INTXBD) == 0
                set(0,'CurrentFigure',fg);
                temp = temp + 1;
                timerecorded(1,temp) = N;
                ZBEDDraw(:,:,temp) = zcont;
                
              % ========= calculate apparent dip =============
              zDiff = diff(z,1,2); % difference between colomn
              xDiff = diff(XSurf,1,2);
              yDiff = diff(YSurf,1,2);

              if min(min(xDiff)) == 0
                  xSlope = yDiff;
              else
                  xSlope = xDiff;
              end
                
              % zSlopeInDgree(0,180] degree, to up 0, horizontal 90, to bottom 180 
              % apparentDip [0,90] degree
              zSlopeInDgree = NaN(size(z,1), size(z,2));
              apparentDip = NaN(size(z,1), size(z,2));
              trueDip = NaN(size(z,1), size(z,2));
              apparentDipSign = NaN(size(z,1), size(z,2));
              for p = 1 : size(zDiff,2)
                  for q = 1 : size(zDiff,1)
                     zSlopethis = vectorInCompassDir(xDiff(q,p),zDiff(q,p)); 
                     zSlopeInDgree(q,p+1) = zSlopethis;
                     apparentDip(q,p+1) = abs(zSlopethis - 90);
                     apparentDipSign(q,p+1) = zSlopethis - 90; % negative: stoss slope; positive: lee slope
                  end
              end
              zSlopeInDgree(:,1) = zSlopeInDgree(:,2);
              apparentDip(:,1) = apparentDip(:,2);
              apparentDipSign(:,1) = apparentDipSign(:,2);
              
              if temp == 1
                  apparentDipDraw(:,:,1) = apparentDip;
                  apparentDipSignDraw(:,:,1) = apparentDipSign;
              else
                  apparentDipDraw(:,:,temp) = apparentDip;
                  apparentDipSignDraw(:,:,temp) = apparentDipSign;
                  
                  tempz = z;
                  [r,c] = find(tempz > zcont); % use the apparent dip from zcont
                  
                  for i = 1 : numel(r)
                      rowthis = r(i);
                      colthis = c(i);

                      apparentDipDraw(rowthis,colthis,temp) = apparentDipDraw(rowthis,colthis,temp-1);
                      apparentDipSignDraw(rowthis,colthis,temp) = apparentDipSignDraw(rowthis,colthis,temp-1);
                      ZErosion(rowthis,colthis,temp) = 1;
                  end
                  
              end
                

                set(0,'CurrentFigure',fg);
                surf(x,y,zcont,'LineStyle','none');
                hold on;
          end

          if max(ZBED) <= -30; MORVRT = false; end
    end

  end  % End of loop to finish each image
  ZBEDDraw(:,:,temp+1:end) = [];
  apparentDipDraw(:,:,temp+1:end) = [];
  apparentDipSignDraw(:,:,temp+1:end) = [];

  % SAVE THE MATRIX FOR HORIZONTAL SLICES
  cd(outputFolder_simulation);
  save('TopoINfO','XSurf','YSurf','ZBEDDraw','apparentDipDraw','apparentDipSignDraw','ZErosion');
  
  
  ZBED_xz = permute(ZBEDDraw,[2 3 1]); % X-by-Z array
  ZBED_yz = permute(ZBEDDraw,[1 3 2]); % Y-by-Z array
  
%   figure;
%   plot(x(1,:),ZBED_xz(:,:,1));
%   
%   figure;
%   plot(y(:,1),ZBED_yz(:,:,1));
  

%   % get cross sections right, from closest
%   figure;
%   for i = 1 : size(ZBED_xz,3)   % i is the specified y value
%       ZBED_xz_this = ZBED_xz(:,:,i);
%       plot3(x(i,:),y(i,:),ZBED_xz_this,'Color', [0.1 .1 .1])     
%       hold on;
%       % drawnow;
%   end
%   hold off;
%   
%   
%   % get cross sections left
%   figure;
%   for j = 1 : size(ZBED_yz,3) % j is specified x direciton
%       
%       ZBED_yz_this = ZBED_yz(:,:,j);
%       plot3(y(j,:),x(j,:),ZBED_yz_this,'Color', [0.1 .1 .1])  
%       hold on;
%   end




%   ImageName = sprintf('Topo_%s_%03.0f', FILENM, slicei);
%   cd(outputFolder_simulation);
%   print('-dtiff','-r250',ImageName)
  % print('-depsc2', '-loose', ImageName);


end

% ------ check data -------
% TIME_list(isnan(TIME_list)) = []; TIME_list(1:2,:) = [];
% DPOSIT_List(isnan(DPOSIT_List)) = []; DPOSIT_List(1:2,:) = [];
% dispfd_List(isnan(dispfd_List)) = []; dispfd_List(1:2,:) = [];
% sizef_List(isnan(sizef_List)) = []; sizef_List(1:2,:) = [];
% figure;
% plot(TIME_list, DPOSIT_List);
% hold on;
% plot(TIME_list, dispfd_List);
% 
% figure;
% plot(TIME_list, sizef_List);


