% calcuate true dip for each cross sections
% clear;
clf
close all;
drawlines = 0;
NoAxisImages = 1;


cd(outputFolder_simulation);
load('TopoINfO.mat');
load('ContourLineINfO.mat');

ZBED_xz = permute(ZBEDDraw,[2 3 1]); % X-by-Z array
ZBED_yz = permute(ZBEDDraw,[1 3 2]); % Y-by-Z array

apparentDipDraw_xz = permute(apparentDipDraw,[2 3 1]);
apparentDipSignDraw_xz = permute(apparentDipSignDraw,[2 3 1]);
ZErosion_xz = permute(ZErosion,[2 3 1]);


% % get 3D surface 
%   figure;
%   for q = 1 : size(ZBEDDraw,3)   % i is the specified y value
%       surf(XSurf,YSurf,ZBEDDraw(:,:,q),'LineStyle','none'); 
%       hold on;
%   end
%   hold off;
  
  
%   % get cross sections left
%   for j = 1 : size(ZBED_yz,3) % j is specified x direciton
%       
%       ZBED_yz_this = ZBED_yz(:,:,j);
%       plot3(y(j,:),x(j,:),ZBED_yz_this,'Color', [0.1 .1 .1])  
% 
%   end

% get cross sections right, from closest
zSlope_allSections = NaN(size(ZBED_xz,1),size(ZBED_xz,2),size(ZBED_xz,3));

for i = 1 : size(ZBED_xz,3)   % i is the specified y value
  
%   set(0,'CurrentFigure',fg);
%   clf;
%   set(gcf,'color','w');  




  
  
  zthis = ZBED_xz(:,:,i);
  xthis = XSurf(i,:);
  ythis = YSurf(i,:);
  ZErosionthis = ZErosion_xz(:,:,i);
  
%   apparentDip = apparentDipDraw_xz(:,:,i);
%   apparentDipSign = apparentDipSignDraw_xz(:,:,i);
  
  % cross section in 3D at ythis
  % figure; plot3(xthis,ythis,zthis,'Color', [0.870588 0.721569 0.529412]); % zthis each coloum is a line
  
  
  % ---- for debug, cross section in 2D ----

  close all;
  f = figure('visible','on');


% figure;

%   for ii = 1 : size(ZBEDDraw,3)
%       plot(xthis',zthis(:,ii),'Color', [0.870588 0.721569 0.529412],'Linewidth',0.5);
%       hold on;
%   end
  % ----------------


  xmatrix = repmat(xthis', 1, size(zthis,2));
  % each coloum each line

  linenumber = size(zthis,2);
  pointsNumberEachLine = size(xthis,2);
  dx = zeros(pointsNumberEachLine,linenumber); 
  dy = zeros(pointsNumberEachLine,linenumber); 
  NVa = zeros(pointsNumberEachLine,linenumber); 
  NVb = zeros(pointsNumberEachLine,linenumber); 


  dx(1,:) = (xmatrix(2,:) - xmatrix(1,:))/1;
  dy(1,:) = (zthis(2,:) - zthis(1,:))/1;
  dx(2:pointsNumberEachLine-1,:) = (xmatrix(3:pointsNumberEachLine,:)- xmatrix(1:pointsNumberEachLine-2,:))/2;
  dy(2:pointsNumberEachLine-1,:) = (zthis(3:pointsNumberEachLine,:)- zthis(1:pointsNumberEachLine-2,:))/2;
  dx(pointsNumberEachLine,:) = (xmatrix(pointsNumberEachLine,:) - xmatrix(pointsNumberEachLine-1,:))/1;
  dy(pointsNumberEachLine,:) = (zthis(pointsNumberEachLine,:) - zthis(pointsNumberEachLine-1,:))/1;


  % zSlopeInDgree(0,180] degree, to up 0, horizontal 90, to bottom 180 
  % apparentDip [0,90] degree
  zSlopeInDgree = NaN(pointsNumberEachLine,linenumber);
  apparentDip = NaN(pointsNumberEachLine,linenumber);
  trueDip = NaN(pointsNumberEachLine,linenumber);
  apparentDipSign = NaN(pointsNumberEachLine,linenumber);
  for p = 1 : linenumber
      for q = 1 : pointsNumberEachLine    
         zSlopethis = vectorInCompassDir(dx(q,p),dy(q,p)); 
         zSlopeInDgree(q,p) = zSlopethis;
         apparentDip(q,p) = abs(zSlopethis - 90);
         apparentDipSign(q,p) = zSlopethis - 90; % negative: stoss slope; positive: lee slope
      end
  end
  % zSlopeInDgree(1,:) = zSlopeInDgree(2,:);
  % apparentDip(1,:) = apparentDip(2,:);
  % apparentDipSign(1,:) = apparentDipSign(2,:);
  
  for j = 1 : size(xthis,2)
     xP = xthis(1,j);
     yP = ythis(1,j);
     zC = zthis(j,:);  % same x and y
     
     for m = 1: size(zthis,2) 
         thisP = zC(1,m);
         absDis = abs(thisP - zrefList); % zrefList, add density
         [M,I] = min(absDis);
         thisHorContLines = ContourLinesList{I}; % ContourLinesList, add density
         
%          % --- draw horizontal section to check ---
%          hold on;
%          for checki = 1 : max(thisHorContLines(:,2))
%              [rcont,~] = find(thisHorContLines(:,2)==checki);
%              thisline = thisHorContLines(rcont,:);
%              plot3(thisline(:,3),thisline(:,4),thisline(:,1));
%              hold on;
%          end
%          plot3(xP,yP,thisP,'ro');
%          % -----------------------------------------
         
         disOnCont = sqrt((thisHorContLines(:,3) - xP).^2 + (thisHorContLines(:,4) - yP).^2);
         [M2,I2] = min(disOnCont);
         thislineNumber = thisHorContLines(I2, 2);
         [rthisline,~] = find(thisHorContLines(:,2) == thislineNumber);
         thisline = thisHorContLines(rthisline,:);
         x = thisline(:,3);
         y = thisline(:,4);
         disP = sqrt((x-xP).^2 + (y-yP).^2);
         [MP,IP] = min(disP);
         if IP > 1 && IP < numel(disP)
             thisP_dir = vectorInCompassDir(x(IP-1)-x(IP+1), y(IP-1)-y(IP+1));
         end
         if IP == 1
             thisP_dir = vectorInCompassDir(x(1)-x(2), y(1)-y(2));
         end  
         if IP == numel(disP)
             thisP_dir = vectorInCompassDir(x(end-1)-x(end), y(end-1)-y(end));
         end
         if thisP_dir <= 180
             strikeDir = thisP_dir;
         else
             strikeDir = thisP_dir - 180;
         end
         
         crossSectionDir = 90;
         strike2apparentDip = abs(strikeDir - crossSectionDir);
         if strike2apparentDip ~=0        
            trueDipthis = atand(tand(apparentDip(j,m))/sind(strike2apparentDip));
            trueDip(j,m) = trueDipthis;
         end

     end
  end
  
  [r,c] = find(apparentDipSign < 0);
  trueDipSign = trueDip;
  for i2 = 1 : numel(r)
      trueDipSign(r(i2),c(i2)) = (-1) * trueDipSign(r(i2),c(i2));
  end
  

  % ============= draw facies ================= 
  % DipDirection
  zflat = min(zthis);
  ZBEDthisFaciesProperties = NaN(size(zthis,1), size(zthis,2));
  for p = 1 : size(zthis,2)
      trueDipSignthis = trueDipSign(:,p);
      zthistime = zthis(:,p);
      zflatthis = zflat(1,p);
      
%       [r0,~] = find(trueDipSignthis == 0 & zthistime == zflatthis);
%       [r1,~] = find(trueDipSignthis <= 8 & zthistime > zflatthis);
%       [r2,~] = find(trueDipSignthis > 8 & trueDipSignthis <= 16);
%       [r3,~] = find(trueDipSignthis > 16 & trueDipSignthis <= 25);
%       [r4,~] = find(trueDipSignthis > 25);

%       dipDefine = [8,16,25]; % boundryies between faces

      [r0,~] = find(trueDipSignthis == 0 & zthistime == zflatthis);
      [r1,~] = find(trueDipSignthis <= dipDefine(1,1) & zthistime > zflatthis);
      [r2,~] = find(trueDipSignthis > dipDefine(1,1) & trueDipSignthis <= dipDefine(1,2));
      [r3,~] = find(trueDipSignthis > dipDefine(1,2) & trueDipSignthis <= dipDefine(1,3));
      [r4,~] = find(trueDipSignthis > dipDefine(1,3));
      
      ZBEDthisFaciesProperties(r0,p) = 0;
      ZBEDthisFaciesProperties(r1,p) = 1;
      ZBEDthisFaciesProperties(r2,p) = 2;
      ZBEDthisFaciesProperties(r3,p) = 3;
      ZBEDthisFaciesProperties(r4,p) = 4;
      
      
  end



%%
  
  for qi = 1 : size(zthis,2)-1
      z1 = zthis(:,qi);
      z2 = zthis(:,qi+1);
      zdiff = z1 - z2;
      
      % find the start and end index of zdiff = 0
      % ii = [0, diff(zdiff(:)')==0,0];

      ii = [0, zdiff(:)'==0,0];
      istart = strfind(ii,[0 1]); % start index
      iend = strfind(ii,[1 0]); % end index
      

      for i3 = 1 : numel(istart)
         z1(istart(i3)+1:iend(i3)-2,:) = NaN; 
         z2(istart(i3)+1:iend(i3)-2,:) = NaN; 
      end
      if ~isempty(istart) && istart(1) == 1
         z1(1,1) = NaN; 
         z2(1,1) = NaN; 
      end
      if ~isempty(istart) && iend(end) == size(z1,1)
         z1(end,1) = NaN; 
         z2(end,1) = NaN; 
      end
      % --------------------

     
      z1facies = ZBEDthisFaciesProperties(:,qi);
      z2facies = ZBEDthisFaciesProperties(:,qi+1);
      z1faciesDipSign = trueDipSign(:,qi);
      z2faciesDipSign = trueDipSign(:,qi+1);
      
      
      z1 = z1';
      z2 = z2';
      z1facies = z1facies';
      z2facies = z2facies';
      z1faciesDipSign = z1faciesDipSign';
      z2faciesDipSign = z2faciesDipSign';
      
      
      index1 = all(~isnan(z1),1);
      idx1 = [strfind([~index1(1),index1],[0 1]); strfind([index1, ~index1(end)],[1 0])];
      z1_parts = mat2cell(z1(:,index1),size(z1,1),diff(idx1)+1);
      z1facies_parts = mat2cell(z1facies(:,index1),size(z1,1),diff(idx1)+1);
      x1_CrossSection_parts = mat2cell(xthis(:,index1),size(z1,1),diff(idx1)+1);
      z1faciesDipSign_parts = mat2cell(z1faciesDipSign(:,index1),size(z1,1),diff(idx1)+1);

      index2 = all(~isnan(z2),1);
      idx2 = [strfind([~index2(1),index2],[0 1]); strfind([index2, ~index2(end)],[1 0])];
      z2_parts = mat2cell(z2(:,index2),size(z2,1),diff(idx2)+1);
      z2facies_parts = mat2cell(z2facies(:,index2),size(z2,1),diff(idx2)+1);
      x2_CrossSection_parts = mat2cell(xthis(:,index2),size(z2,1),diff(idx2)+1);
      z2faciesDipSign_parts = mat2cell(z2faciesDipSign(:,index2),size(z2,1),diff(idx2)+1);

      partsN = min(size(z1_parts,2),size(z2_parts,2));
      
      
      
      for ip = 1 : partsN 
          
          if ip == 2
             a = 1; 
          end
          
          z1this = z1_parts{ip}; x1this = x1_CrossSection_parts{ip};
          z2this = z2_parts{ip}; x2this = x2_CrossSection_parts{ip};
          z1dipSign = z1faciesDipSign_parts{ip};
          z2dipSign = z2faciesDipSign_parts{ip};
          
          z1faciesthis = z1facies_parts{ip};       
          z2faciesthis = z2facies_parts{ip}; 
          
          % debug
%           figure;
%           plot(x1this,z1this,'o');
% 
%           A = [x1this',z1this'];
%           k = LineCurvature2D(A);
%           figure;
%           plot(k);

          
          
          hold on;
          
          if numel(x1this)>2
          B1 = smoothdata([x1this',z1this'],"movmean",3); 
          B2 = smoothdata([x2this',z2this'],"movmean",3);  
          TF1 = islocalmax(B1(:,2)','SamplePoints',B1(:,1)');
          TF3 = islocalmin(B1(:,2)','SamplePoints',B1(:,1)');
          TF2 = islocalmax(B2(:,2)','SamplePoints',B2(:,1)');
          TF4 = islocalmin(B2(:,2)','SamplePoints',B2(:,1)');   
          else



          TF1 = islocalmax(z1this,'SamplePoints',x1this);
          TF3 = islocalmin(z1this,'SamplePoints',x1this);

          % plot(x1this,z1this,'r-', x1this(TF1),z1this(TF1),'r*');
          
          TF2 = islocalmax(z2this,'SamplePoints',x2this);
          TF4 = islocalmin(z2this,'SamplePoints',x2this);  
          % plot(x2this,z2this,'b-', x2this(TF2),z2this(TF2),'b*');
          end



          [~,TF1_index] = find(TF1 == 1);
          [~,TF2_index] = find(TF2 == 1);
          [~,TF3_index] = find(TF3 == 1);
          [~,TF4_index] = find(TF4 == 1);

          

          % bottom

          BTB = 0;
          if numel(TF3_index) > 1 && numel(TF4_index) > 1    % TF3_index(local min, x1this,z1this); TF4_index(local min, x2this,z2this)
               % [~,c1] = find(TF3_index > 3);
               % [~,c2] = find(TF4_index > 3);
               
               %if numel(c1) > 1 && numel(c2) > 1
                    BTB = 1;


                    if numel(TF1_index) == 1
                        [~,indextemp] = min(abs(TF2_index - TF1_index));
                        TF2_index = TF2_index(1,indextemp);

                    elseif numel(TF2_index) == 1
                        [~,indextemp] = min(abs(TF2_index - TF1_index));
                        TF1_index = TF1_index(1,indextemp);

                    else

                    % TF1_index = TF1_index(1); % TF1_index(local max, x1this,z1this); 
                    % TF2_index = TF2_index(1); % TF2_index(local max, x2this,z2this); 

                    colmid = 0.5 * numel(z1this);
                    k1 = dsearchn(TF1_index', colmid);
                    TF1_index = TF1_index(k1);

                    k2 = dsearchn(TF2_index', colmid);
                    TF2_index = TF2_index(k2);


                    end


                    % TF3_index, TF4_index need only 1


                    
               %end

          end


          % only 1 bottom
          if BTB == 0 
              if numel(TF3_index) >= 1
                    [~,c] = find(TF3_index > 3);
     
                    if numel(c) > 1
                        if TF3_index(1,c(1)) > numel(z1this) - TF3_index(1,c(end))
                            TF3_index = TF3_index(1,c(1));
                        else
                            TF3_index = TF3_index(1,c(end));
                        end
                    
                    elseif numel(c) == 1
                        TF3_index = TF3_index(1,c);
    
                    else
                        TF3_index = [];
                    end
              end
    
              if numel(TF4_index) >= 1
                    [~,c] = find(TF4_index > 3);
                    
                    if numel(c) > 1
                        if TF4_index(1,c(1)) > numel(z2this) - TF4_index(1,c(end))
                            TF4_index = TF4_index(1,c(1));
                        else
                            TF4_index = TF4_index(1,c(end));
                        end
    
                    elseif numel(c) == 1
%                         if abs(TF4_index(1,c) - TF3_index) > (max(z1this) - min(z1this))*0.3
%                             TF4_index = NaN;
%                         else
%                             TF4_index = TF4_index(1,c);
%                         end

                        TF4_index = TF4_index(1,c);

                        
    
                    else
                        TF4_index = [];
                    end
              end
          end


          



%           % ----------- debug, check point --------
%           % cross section in 3D at ythis
%           plot3(xthis,ythis,zthis,'Color', [0.870588 0.721569 0.529412]); % zthis each coloum is a line
% 
%           plot3(x1this,ythis(1:size(z1this,2)),z1this,'Color', [1 0 0]); 
%           if sum(TF1)>0 
%               plot3(x1this(TF1),ythis(1:size(TF1,2)), z1this(TF1),'Color', [1 0 0],'Marker','o') 
%           end
%           plot3(x2this,ythis(1:size(x2this,2)),z2this,'Color', [0 0 1]); % zthis each coloum is a line
%           if sum(TF2)>0
%               plot3(x2this(TF2),ythis(1:size(TF2,2)), z2this(TF2),'Color', [0 0 1],'Marker','*');
%           end
% 
%           checkpointNumber = 36;
%           if checkpointNumber<= size(z1this,2)
%              zplot = z1this(checkpointNumber);       
%              absDis = abs(zplot - zrefList); % zrefList, add density
%              [M,I] = min(absDis);
%              thisHorContLines = ContourLinesList{I}; % ContourLinesList, add density
%          
%              % --- draw horizontal section to check ---
%              hold on;
%              for checki = 1 : max(thisHorContLines(:,2))
%                  [rcont,~] = find(thisHorContLines(:,2)==checki);
%                  thisline = thisHorContLines(rcont,:);
%                  plot3(thisline(:,3),thisline(:,4),thisline(:,1), 'Color',[0.5 0.5 0.5]);
%                  hold on;
%              end
%              plot3(x1this(checkpointNumber),ythis(1),z1this(checkpointNumber),'m^');
%           end
%           % -----------------------------------------
         
               
            

          % ---------- deal with nan of facies, use the point before ----
          [~,colnan1] = find(isnan(z1faciesthis));
          if ~isempty(colnan1)
              for i4 = 1 : numel(colnan1)
                  if colnan1(i4)>1
                     temp = 1;
                     while colnan1(i4)-temp>0 && isnan(z1faciesthis(colnan1(i4)-temp))
                        temp = temp + 1;
                     end
                     if colnan1(i4)-temp>0
                        z1faciesthis(colnan1(i4)) = z1faciesthis(colnan1(i4)-temp);
                     end
                  end
              end
          end
          
          [~,colnan3] = find(isnan(z1faciesthis));
          if ~isempty(colnan3)
              for i4 = 1 : numel(colnan3)
                 temp = 1;
                 while colnan3(i4)+temp <= numel(z1faciesthis) && isnan(z1faciesthis(colnan3(i4)+temp))
                    temp = temp + 1;
                 end
                 if colnan3(i4)+temp <= numel(z1faciesthis)
                    z1faciesthis(colnan3(i4)) = z1faciesthis(colnan3(i4)+temp);
                 end
              end
          end
          
          [~,colnan2] = find(isnan(z2faciesthis));
          if ~isempty(colnan2)
              for i5 = 1 : numel(colnan2)
                  if colnan2(i5)>1
                     temp = 1;
                     while colnan2(i5)-temp>0 && isnan(z2faciesthis(colnan2(i5)-temp))
                        temp = temp + 1;
                     end
                     if colnan2(i5)-temp>0
                         z2faciesthis(colnan2(i5)) = z2faciesthis(colnan2(i5)-temp);
                     end
                  end
              end
          end
          
           [~,colnan4] = find(isnan(z2faciesthis));
          if ~isempty(colnan4)
              for i5 = 1 : numel(colnan4)
                 temp = 1;
                 while colnan4(i5)+temp <= numel(z2faciesthis) && isnan(z2faciesthis(colnan4(i5)+temp))
                    temp = temp + 1;
                 end
                 if colnan4(i5)+temp <= numel(z2faciesthis)
                    z2faciesthis(colnan4(i5)) = z2faciesthis(colnan4(i5)+temp);
                 end
              end
          end
          
          
          z1faciesthis(isnan(z1faciesthis)) = [];
          z2faciesthis(isnan(z2faciesthis)) = [];
          if isempty(z1faciesthis) || isempty(z2faciesthis)
              continue;
          end
          
          
          
          
          z1faciesNumber = z1faciesthis;
          z1faciesNumber(diff(z1faciesNumber)==0) = [];
          z2faciesNumber = z2faciesthis;
          z2faciesNumber(diff(z2faciesNumber)==0) = [];
          
          z1faciesthis(1,end+1) = 10; % for format reason, any useless number can do
          I_end1 = find(diff(z1faciesthis)~=0);
          if numel(I_end1) > 1
            I_start1 = I_end1(1,1:end-1) + 1;
            I_start1 = [1,I_start1];
          else
            I_start1 = 1;
          end
          z1faciesIndex = [z1faciesNumber',I_start1',I_end1'];
          z1faciesIndex_x = [x1this(1,z1faciesIndex(:,2)');x1this(1,z1faciesIndex(:,3)')];
          z1faciesIndex_z = [z1this(1,z1faciesIndex(:,2)');z1this(1,z1faciesIndex(:,3)')];
          z1faciesIndex_dipSign = [z1dipSign(1,z1faciesIndex(:,2)');z1dipSign(1,z1faciesIndex(:,3)')];
          z1faciesIndexNew = [z1faciesIndex, z1faciesIndex_x',z1faciesIndex_z',z1faciesIndex_dipSign'];
          % z1faciesIndexNew [facies type, point start, point end, x start, x_end, z_start,
          % z_end, dipsign_start, dipsign_end]
          
%           % ---- for debugging ----
%           figure;
%           plot(x1this, z1this,'*-'); % acceration surface line
%           hold on;
%           plot(z1faciesIndexNew(:,4),z1faciesIndexNew(:,6),'ro'); % start point of each facies
%           plot(z1faciesIndexNew(:,5),z1faciesIndexNew(:,7),'b^'); % end point of each facies
% 
%           figure;
%           plot(z1dipSign);
%           % -----------------------
         
          
          if ~isempty(z2faciesthis)
          z2faciesthis(1,end+1) = 10; % for format reason, any useless number can do
          I_end = find(diff(z2faciesthis)~=0);
          if numel(I_end) > 1
            I_start = I_end(1,1:end-1) + 1;
            I_start = [1,I_start];
          else
            I_start = 1;
          end
          z2faciesIndex = [z2faciesNumber',I_start',I_end'];
          z2faciesIndex_x = [x2this(1,z2faciesIndex(:,2)');x2this(1,z2faciesIndex(:,3)')];
          z2faciesIndex_z = [z2this(1,z2faciesIndex(:,2)');z2this(1,z2faciesIndex(:,3)')];
          z2faciesIndex_dipSign = [z2dipSign(1,z2faciesIndex(:,2)');z2dipSign(1,z2faciesIndex(:,3)')];
          z2faciesIndexNew = [z2faciesIndex, z2faciesIndex_x',z2faciesIndex_z',z2faciesIndex_dipSign'];
          else
            z2faciesIndexNew = [];
          end

          % z2faciesIndexNew [facies type, point start, point end, x start, x_end, z_start,
          % z_end, dipsign_start, dipsign_end]
          
          %##########  x1this, x2this, z1this, z2this, z1faciesthis,
          %z2faciesthis, z1faciesIndexNew, z2faciesIndexNew
          
          % <=8 [1]
          % >8 <=16 [2]
          % >16 <=25 [3]
          % >25 [4]

          % switch
          xthistemp1 = x1this;
          zthistemp1 = z1this;
          zfaciesthistemp1 = z1faciesthis;
          zfaciesindextemp1 = z1faciesIndexNew;

          xthistemp2 = x2this;
          zthistemp2 = z2this;
          zfaciesthistemp2 = z2faciesthis;
          zfaciesindextemp2 = z2faciesIndexNew;

          x1this = xthistemp2;
          z1this = zthistemp2;
          z1faciesthis = zfaciesthistemp2;
          z1faciesIndexNew = zfaciesindextemp2;

          x2this = xthistemp1;
          z2this = zthistemp1;
          z2faciesthis = zfaciesthistemp1;
          z2faciesIndexNew = zfaciesindextemp1;


          % max of line
          tf1 = TF1_index;
          tf2 = TF2_index;
          TF1_index = tf2;
          TF2_index = tf1;

          % min of line
          tf3 = TF3_index;
          tf4 = TF4_index;
          TF3_index = tf4;
          TF4_index = tf3;



          % facies type only include one points, delete
          if size(z1faciesIndexNew,1) > 1

              z1faciesIndexNew2 = z1faciesIndexNew;
              for m = size(z1faciesIndexNew2,1):-1:1
                 if z1faciesIndexNew2(m,2) == z1faciesIndexNew2(m,3)
    
                    if m == 1
                        z1faciesIndexNew2(m+1,2) = z1faciesIndexNew2(m,3);
                        z1faciesIndexNew2(m+1,4) = z1faciesIndexNew2(m,5);
                        z1faciesIndexNew2(m+1,6) = z1faciesIndexNew2(m,7);
                        z1faciesIndexNew2(m+1,8) = z1faciesIndexNew2(m,9);
                        z1faciesIndexNew2(m,:) = [];
    
                    elseif m == size(z1faciesIndexNew2,1)
                        z1faciesIndexNew2(m-1,3) = z1faciesIndexNew2(m,2);
                        z1faciesIndexNew2(m-1,5) = z1faciesIndexNew2(m,4);
                        z1faciesIndexNew2(m-1,7) = z1faciesIndexNew2(m,6);
                        z1faciesIndexNew2(m-1,9) = z1faciesIndexNew2(m,8);
                        z1faciesIndexNew2(m,:) = [];
                    else
                        zcheck = z1faciesIndexNew2(m,6);
                        zbefore = z1faciesIndexNew2(m-1,7);
                        zafter = z1faciesIndexNew2(m+1,6);
    
                        if abs(zcheck-zbefore) <= abs(zcheck-zafter) || z1faciesIndexNew2(m+1,1) == 0
                            z1faciesIndexNew2(m-1,3) = z1faciesIndexNew2(m,2);
                            z1faciesIndexNew2(m-1,5) = z1faciesIndexNew2(m,4);
                            z1faciesIndexNew2(m-1,7) = z1faciesIndexNew2(m,6);
                            z1faciesIndexNew2(m-1,9) = z1faciesIndexNew2(m,8);
                            z1faciesIndexNew2(m,:) = [];
                        else
                            z1faciesIndexNew2(m+1,2) = z1faciesIndexNew2(m,3);
                            z1faciesIndexNew2(m+1,4) = z1faciesIndexNew2(m,5);
                            z1faciesIndexNew2(m+1,6) = z1faciesIndexNew2(m,7);
                            z1faciesIndexNew2(m+1,8) = z1faciesIndexNew2(m,9);
                            z1faciesIndexNew2(m,:) = [];
                        end
    
                    end
                 end
    
              end

              z1faciesIndexNew1 = z1faciesIndexNew;

              m = 1;
              while (m < size(z1faciesIndexNew1,1))
              % for m = 1: size(z1faciesIndexNew1,1)
                 if z1faciesIndexNew1(m,2) == z1faciesIndexNew1(m,3)
    
                    if m == 1
                        z1faciesIndexNew1(m+1,2) = z1faciesIndexNew1(m,3);
                        z1faciesIndexNew1(m+1,4) = z1faciesIndexNew1(m,5);
                        z1faciesIndexNew1(m+1,6) = z1faciesIndexNew1(m,7);
                        z1faciesIndexNew1(m+1,8) = z1faciesIndexNew1(m,9);
                        z1faciesIndexNew1(m,:) = [];
                        m = 1;
                       
    
                    elseif m == size(z1faciesIndexNew1,1)
                        z1faciesIndexNew1(m-1,3) = z1faciesIndexNew1(m,2);
                        z1faciesIndexNew1(m-1,5) = z1faciesIndexNew1(m,4);
                        z1faciesIndexNew1(m-1,7) = z1faciesIndexNew1(m,6);
                        z1faciesIndexNew1(m-1,9) = z1faciesIndexNew1(m,8);
                        z1faciesIndexNew1(m,:) = [];
             

                    else
                        zcheck = z1faciesIndexNew1(m,6);
                        zbefore = z1faciesIndexNew1(m-1,7);
                        zafter = z1faciesIndexNew1(m+1,6);
    
                        if abs(zcheck-zbefore) <= abs(zcheck-zafter) || z1faciesIndexNew1(m+1,1) == 0
                            z1faciesIndexNew1(m-1,3) = z1faciesIndexNew1(m,2);
                            z1faciesIndexNew1(m-1,5) = z1faciesIndexNew1(m,4);
                            z1faciesIndexNew1(m-1,7) = z1faciesIndexNew1(m,6);
                            z1faciesIndexNew1(m-1,9) = z1faciesIndexNew1(m,8);
                            z1faciesIndexNew1(m,:) = [];

                        else
                            z1faciesIndexNew1(m+1,2) = z1faciesIndexNew1(m,3);
                            z1faciesIndexNew1(m+1,4) = z1faciesIndexNew1(m,5);
                            z1faciesIndexNew1(m+1,6) = z1faciesIndexNew1(m,7);
                            z1faciesIndexNew1(m+1,8) = z1faciesIndexNew1(m,9);
                            z1faciesIndexNew1(m,:) = [];
                        end
    
                    end
                 end
                 m = m + 1;
    
              end



          z1faciesIndexNew = z1faciesIndexNew1;
          end
          


          % merge the boundary between two facies
          for m1 = 1:size(z1faciesIndexNew,1)-1
              dip1 = z1faciesIndexNew(m1,9);
              dip2 = z1faciesIndexNew(m1+1,8);
              dipdiff1 = abs(dipDefine - abs(dip1));
              dipdiff2 = abs(dipDefine - abs(dip2));
              if min(dipdiff1) <= min(dipdiff2)
                  z1faciesIndexNew(m1+1,2) = z1faciesIndexNew(m1,3);
                  z1faciesIndexNew(m1+1,4) = z1faciesIndexNew(m1,5);
                  z1faciesIndexNew(m1+1,6) = z1faciesIndexNew(m1,7);
                  z1faciesIndexNew(m1+1,8) = z1faciesIndexNew(m1,9);

              else
                  z1faciesIndexNew(m1,3) = z1faciesIndexNew(m1+1,2);
                  z1faciesIndexNew(m1,5) = z1faciesIndexNew(m1+1,4);
                  z1faciesIndexNew(m1,7) = z1faciesIndexNew(m1+1,6);
                  z1faciesIndexNew(m1,9) = z1faciesIndexNew(m1+1,8);



              end
          end


          
          % ------------ MATCHING TOP AND BOTTOM IN ORDER -------

          if BTB == 0
          TBT_1 = NaN(1,4);  
          if ~isempty(TF3_index)
                TBT_1(1,2) = TF3_index; % bottom
          end
          if ~isempty(TF1_index)
              for ti = 1 : numel(TF1_index)
                    if isempty(TF3_index)
                        if TF1_index(1,ti) < 0.5 * numel(x1this)
                            TBT_1(1,1) = TF1_index(1,ti);
                        else
                            TBT_1(1,3) = TF1_index(1,ti);
                        end
                        
                    elseif ~isempty(TF3_index) && TF1_index(1,ti) < TF3_index
                        TBT_1(1,1) = TF1_index(1,ti);
                    elseif ~isempty(TF3_index) && TF1_index(1,ti) > TF3_index
                        TBT_1(1,3) = TF1_index(1,ti);
                    else % TF3_index is empty
                        TBT_1(1,1) = TF1_index(1,ti);
                    end

              end
          end
          % TBT_1(1,4) = numel(x1this);

          TBT_2 = NaN(1,4);  
          if ~isempty(TF4_index)
                TBT_2(1,2) = TF4_index; % top
          end
          if ~isempty(TF2_index)
              for ti = 1 : numel(TF2_index)
                    if isempty(TF4_index)
                        if TF2_index(1,ti) < 0.5 * numel(x2this)
                            TBT_2(1,1) = TF2_index(1,ti);
                        else
                            TBT_2(1,3) = TF2_index(1,ti);
                        end


                        
                    elseif ~isempty(TF4_index) && TF2_index(1,ti) < TF4_index
                        TBT_2(1,1) = TF2_index(1,ti);
                    elseif ~isempty(TF4_index) && TF2_index(1,ti) > TF4_index
                        TBT_2(1,3) = TF2_index(1,ti);
                    else

                    end

              end
          end
          % TBT_2(1,4) = numel(x2this);
          end

          if BTB == 1
          TBT_1 = NaN(1,4);  
          if ~isempty(TF1_index)
                TBT_1(1,2) = TF1_index; % top
          end
          if ~isempty(TF3_index)
              for ti = 1 : numel(TF3_index)
                    if TF3_index(1,ti) < TF1_index
                        TBT_1(1,1) = TF3_index(1,ti);
                    else
                        TBT_1(1,3) = TF3_index(1,ti);
                    end

              end
          end
          % TBT_1(1,4) = numel(x1this);

          TBT_2 = NaN(1,4);  
          if ~isempty(TF2_index)
                TBT_2(1,2) = TF2_index;
          end
          if ~isempty(TF4_index)
              for ti = 1 : numel(TF4_index)
                    if isempty(TF4_index)
                        if TF4_index(1,ti) < 0.5 * numel(x2this)
                            TBT_2(1,1) = TF4_index(1,ti);
                        else
                            TBT_2(1,3) = TF4_index(1,ti);
                        end

                    elseif ~isempty(TF2_index) && TF4_index(1,ti) < TF2_index
                        TBT_2(1,1) = TF4_index(1,ti);
                    elseif ~isempty(TF2_index) && TF4_index(1,ti) > TF2_index
                        TBT_2(1,3) = TF4_index(1,ti);
                    else

                    end

              end
          end
          % TBT_2(1,4) = numel(x2this);

          end


          TBT_1(1,4) = numel(x1this);
          TBT_2(1,4) = numel(x2this);

          % deal with the beginning of curve
          if size(z1faciesIndexNew,1)>1 && z1faciesIndexNew(1,3) <= 3 && z1faciesIndexNew(1,1)~=0
                z1faciesIndexNew(2,2) = 1;
                z1faciesIndexNew(2,4) = z1faciesIndexNew(1,4);
                z1faciesIndexNew(2,6) = z1faciesIndexNew(1,6);
                z1faciesIndexNew(2,8) = z1faciesIndexNew(1,8);
                z1faciesIndexNew(1,:) = [];
          end



          % make sure local maxs and mins are matching
          dis12 = abs(TBT_1 - TBT_2);
          for tempthisi = 1 : 3
              if dis12(tempthisi) >  0.3 * numel(z1this)
                  TBT_1(tempthisi) = NaN;
                  TBT_2(tempthisi) = NaN;
              end

              % big diff of z
              if ~isnan(TBT_1(tempthisi)) && ~isnan(TBT_2(tempthisi)) && abs(z1this(1,TBT_1(tempthisi)) - z2this(1,TBT_2(tempthisi))) > 0.7 * (max(z1this)-min(z1this))
                  TBT_1(tempthisi) = NaN;
                  TBT_2(tempthisi) = NaN;

              end
          end
          % June 2023
          sumstemp = TBT_1 + TBT_2;
          [~,nanc] = find(isnan(sumstemp));
          TBT_1(:,nanc) = NaN;
          TBT_2(:,nanc) = NaN;



          % delete the same facies line next to each other
          faciestemp = find(diff(z1faciesIndexNew(:,1))==0);
          for i4 = numel(faciestemp):-1:1
              colthis = faciestemp(i4);
              z1faciesIndexNew(colthis,3) = z1faciesIndexNew(colthis+1,3);
              z1faciesIndexNew(colthis,5) = z1faciesIndexNew(colthis+1,5);
              z1faciesIndexNew(colthis,7) = z1faciesIndexNew(colthis+1,7);
              z1faciesIndexNew(colthis,9) = z1faciesIndexNew(colthis+1,9);
              z1faciesIndexNew(colthis+1,:) = [];

          end


          zindex = NaN(size(z1faciesIndexNew,1)-1,1);
          for i3 = 1 : size(z1faciesIndexNew,1) - 1
              zcheck = z1faciesIndexNew(i3,7); % first boundary between two facies
              xcheck = z1faciesIndexNew(i3,5);
              colcheck = z1faciesIndexNew(i3,3);
              
              
             
              for ci = 1 : 4
                   if ~isnan(TBT_1(1,ci)) && xcheck <= x1this(1,TBT_1(1,ci)) % change TBT_1 to TBT_2
                   % if ~isnan(TBT_2(1,ci)) && xcheck <= x2this(1,TBT_2(1,ci)) % change TBT_1 to TBT_2
                       
                       if ~isnan(TBT_2(1,ci))
                           if ci == 1
                               zpart = z2this(1:TBT_2(1,ci)); xpart = x2this(1:TBT_2(1,ci));
                               zIsupp = 0;
                           elseif TBT_2(1,ci) == 1
                               zpart = z2this(1:TBT_2(1,ci)); xpart = x2this(1:TBT_2(1,ci));
                               zIsupp = 0;

                           elseif TBT_2(1,ci) > 1
                               thisci = ci - 1;
                               while isnan(TBT_2(1,thisci)) && thisci > 1
                                   thisci = thisci - 1;
                                   if thisci == 1
                                        break;
                                   end
                               end

                               if ~isnan(TBT_2(1,thisci))
                                   zpart = z2this(TBT_2(1,thisci):TBT_2(1,ci)); xpart = x2this(TBT_2(1,thisci):TBT_2(1,ci));
                                   zIsupp = TBT_2(1,thisci)-1;
                               else % nan for thisci=1
                                   zpart = z2this(1:TBT_2(1,ci)); xpart = x2this(1:TBT_2(1,ci));
                                   zIsupp = 0;
                               end
               

                           else
                           end

                       else % isnan(TBT_2(1,ci))
                           
                           % find end boundary
                           thisci2 = ci+1;
                           while isnan(TBT_2(1,thisci2)) && thisci2 < 4
                                   thisci2 = thisci2 + 1;
                                   if thisci2 == 4
                                        break;
                                   end
                           end

                           % find start boundary
                           thisci3 = thisci2 - 1;
                           while isnan(TBT_2(1,thisci3)) && thisci3 > 1
                               thisci3 = thisci3 - 1;
                               if thisci3 == 1
                                    break;
                               end
                           end
                           if ~isnan(TBT_2(1,thisci3))
                               zpart = z2this(TBT_2(1,thisci3):TBT_2(1,thisci2)); xpart = x2this(TBT_2(1,thisci3):TBT_2(1,thisci2));
                               zIsupp = TBT_2(1,thisci3)-1;
                           else % nan for thisci3=1
                               zpart = z2this(1:TBT_2(1,thisci2)); xpart = x2this(1:TBT_2(1,thisci2));
                               zIsupp = 0;
                           end


                       end


                   break;
                   end
                   
              end


              % 14 June 2023, catch no max and min 
              if isnan(TBT_1(1,1)) && isnan(TBT_1(1,2)) && isnan(TBT_1(1,3))
                  zpart = z2this;
                  xpart = x2this;
                  zIsupp = 0;
              end
     

              dis_zpart = abs(zpart - zcheck);
              [~,zI1] = find(dis_zpart == min(dis_zpart));


              % 10 June 2023
              if numel(zI1) > 1
                  dis2_zpart = (zpart(1,zI1) - zcheck).^2;
                  dis2_xpart = (xpart(1,zI1) - xcheck).^2;
                  dis2_all = sqrt(dis2_zpart + dis2_xpart);
                  [~,disI] = find(dis2_all == min(dis2_all));
                  zI1this = zI1(1,disI);
                  
              else
                  zI1this = zI1;

              end
    
              if numel(zI1this) > 1
                  discol = abs(zI1this - colcheck);
                  [~,disIcol] = find(discol == min(discol));
                  zI1this = zI1this(1,disIcol(end));
              end
    
              
              zindex(i3,1) = zI1this + zIsupp; % point number on line


              
              if i3+1<=size(z1faciesIndexNew,1) &&  z1faciesIndexNew(i3+1,1)==0 % flat interdune surface

                    [r0,~] = find(z2faciesIndexNew(:,1)==0);

                    if numel(r0)==1
                        zindex(i3,1) = z2faciesIndexNew(r0,2);
                    elseif isempty(r0)
                        % zindex(i3,1) = z2faciesIndexNew(end,3);
                    else
                        % do nothing, keep the value from above
                    end

              

              end
          end


          zindex = [1; zindex; numel(z1this)];

          
          % check order small-large
          zindex = sortrows(zindex);
          
          % z1faciesIndexNew, zindex
          % z2faciesIndexNew [facies type, point start, point end, x start, x_end, z_start,
          % z_end, dipsign_start, dipsign_end]


          

          gruopn = size(z1faciesIndexNew,1);

          kall = LineCurvature2D([x1this',z1this']); 
          k2all = LineCurvature2D([x2this',z2this']); 
          for ji = 1 : gruopn
              
              zleft = z1this(z1faciesIndexNew(ji,2):z1faciesIndexNew(ji,3));
              xleft = x1this(z1faciesIndexNew(ji,2):z1faciesIndexNew(ji,3));
              dipleft = z1dipSign(z1faciesIndexNew(ji,2):z1faciesIndexNew(ji,3));

              zright = z2this(zindex(ji,1):zindex(ji+1,1));
              xright = x2this(zindex(ji,1):zindex(ji+1,1));
              dipright = z2dipSign(zindex(ji,1):zindex(ji+1,1));

              zthisPolygon = [zleft, fliplr(zright),zleft(1)];
              xthisPolygon = [xleft, fliplr(xright),xleft(1)];



%               A = [xleft', zleft']; B = [xright',zright'];
%               k = LineCurvature2D(A); k2 = LineCurvature2D(B);

              k = kall(z1faciesIndexNew(ji,2):z1faciesIndexNew(ji,3));
              k2 = k2all(zindex(ji,1):zindex(ji+1,1));




              InflectionPointsXYCol_A = InflectionPoints(xleft, zleft, size(xleft,2)); % [x,y,col]  
              InfelctionPointNumberA = size(InflectionPointsXYCol_A,1);

              InflectionPointsXYCol_B = InflectionPoints(xright, zright, size(xright,2)); % [x,y,col]  
              InfelctionPointNumberB = size(InflectionPointsXYCol_B,1);



              locAreaTop = 0;
              if InfelctionPointNumberA == 0
                  [r1,~] = find(k<0);
                  [r2,~] = find(k>0);
                  if numel(r1)>=numel(r2) % top
                      locAreaTop = 1;
                  else %  bottom
                      locAreaTop = 0;
                  end


              elseif InfelctionPointNumberA == 1 && InfelctionPointNumberB == 0

                  if ~isnan(k2)
                    [r3,~] = find(k2<0);
                    [r4,~] = find(k2>0);
                  else
                    [r3,~] = find(k<0);
                    [r4,~] = find(k>0);

                  end
                  if numel(r3)>=numel(r4) % top
                      locAreaTop = 1;
                  else %  bottom
                      locAreaTop = 0;
                  end



              elseif InfelctionPointNumberA == 1 && InfelctionPointNumberB >= 1

                  locAreaTop = 2;
                  inflectionPointColA = InflectionPointsXYCol_A(1,3);
                  a = abs(InflectionPointsXYCol_B(:,2) - InflectionPointsXYCol_A(1,2));
                  [r,c] = find(a == min(a(:)));
                  inflectionPointColB = InflectionPointsXYCol_B(r,3);
                  if k(1,1) < 0
                      topz_L = zleft(1:inflectionPointColA);
                      topx_L = xleft(1:inflectionPointColA);
                      bottomz_L = zleft(inflectionPointColA:end);
                      bottomx_L = xleft(inflectionPointColA:end);
    
                      topz_R = zright(1:inflectionPointColB);
                      topx_R = xright(1:inflectionPointColB);
                      bottomz_R = zright(inflectionPointColB:end);
                      bottomx_R = xright(inflectionPointColB:end);
    

                  else
                      topz_L = zleft(inflectionPointColA:end);
                      topx_L = xleft(inflectionPointColA:end);
                      bottomz_L = zleft(1:inflectionPointColA);
                      bottomx_L = xleft(1:inflectionPointColA);
    
                      topz_R = zright(inflectionPointColB:end);
                      topx_R = xright(inflectionPointColB:end);
                      bottomz_R = zright(1:inflectionPointColB);
                      bottomx_R = xright(1:inflectionPointColB);
                  end
                  zthisPolygon_top = [topz_L, fliplr(topz_R),topz_L(1)];
                  xthisPolygon_top = [topx_L, fliplr(topx_R),topx_L(1)];

                  zthisPolygon_bottom = [bottomz_L, fliplr(bottomz_R),bottomz_L(1)];
                  xthisPolygon_bottom = [bottomx_L, fliplr(bottomx_R),bottomx_L(1)];


              elseif InfelctionPointNumberA > 1
                  if sum(k) < 0
                      locAreaTop = 1;
                  else
                      locAreaTop = 0;
                  end

              else

              end

              
           
              
              
              

                  faciesname = z1faciesIndexNew(ji,1);
                  [r4,c4] = find(z1faciesIndexNew(:,1) == 4);

                  if ~isempty(r4)
                      facies4 = z1faciesIndexNew(r4,:);
                      facies4pns = facies4(:,3) - facies4(:,2);
                      [V4pns,I4pns] = max(facies4pns);
                      rowNumofFacies4 = r4(I4pns);

                  end




                  switch(faciesname)
                      case 0
                          fill(xthisPolygon,zthisPolygon,[37/256,111/256, 230/256],'EdgeColor','none','FaceAlpha',1); 
                      case 4  
                          fill(xthisPolygon,zthisPolygon,[255/256,237/256, 0/256],'EdgeColor','none','FaceAlpha',1); 
                      case 3
                          if locAreaTop == 0
                              fill(xthisPolygon,zthisPolygon,[255/256,191/256, 0/256],'EdgeColor','none','FaceAlpha',1); 
                          elseif locAreaTop == 1
                              fill(xthisPolygon,zthisPolygon,[247/256,171/256, 110/256],'EdgeColor','none','FaceAlpha',1); 
                             
                          elseif locAreaTop == 2
                              fill(xthisPolygon_bottom,zthisPolygon_bottom,[255/256,191/256, 0/256],'EdgeColor','none','FaceAlpha',1); 
                              fill(xthisPolygon_top,zthisPolygon_top,[247/256,171/256, 110/256],'EdgeColor','none','FaceAlpha',1); 
                          else

                          end

                      case 2

                          % fill(xthisPolygon,zthisPolygon,[140/256,115/256, 114/256],'EdgeColor','none','FaceAlpha',1); 
                          % locAreaTop = 0;
                          if locAreaTop == 0
                              fill(xthisPolygon,zthisPolygon,[140/256,115/256, 114/256],'EdgeColor','none','FaceAlpha',1); 
                          elseif locAreaTop == 1
                              fill(xthisPolygon,zthisPolygon,[247/256,171/256, 110/256],'EdgeColor','none','FaceAlpha',1); 
                             
                          elseif locAreaTop == 2
                              fill(xthisPolygon_bottom,zthisPolygon_bottom,[140/256,115/256, 114/256],'EdgeColor','none','FaceAlpha',1); 
                              fill(xthisPolygon_top,zthisPolygon_top,[247/256,171/256, 110/256],'EdgeColor','none','FaceAlpha',1); 
                          else

                          end


                      case 1
                          

                          
                               % fill(xthisPolygon,zthisPolygon,[96/256,64/256, 44/256],'EdgeColor','none','FaceAlpha',1); 
                          % locAreaTop = 0;
                          if locAreaTop == 0
                              fill(xthisPolygon,zthisPolygon,[96/256,64/256, 44/256],'EdgeColor','none','FaceAlpha',1); 
                          elseif locAreaTop == 1
                              fill(xthisPolygon,zthisPolygon,[247/256,171/256, 110/256],'EdgeColor','none','FaceAlpha',1); 
                             
                          elseif locAreaTop == 2
                              fill(xthisPolygon_bottom,zthisPolygon_bottom,[96/256,64/256, 44/256],'EdgeColor','none','FaceAlpha',1); 
                              fill(xthisPolygon_top,zthisPolygon_top,[247/256,171/256, 110/256],'EdgeColor','none','FaceAlpha',1); 
                          else

                          end

                  end

              

          end
          hold on;
          
          
      end
      
  end
  
  if drawlines == 1
      for ii = 1 : size(ZBEDDraw,3)
          plot(xthis',zthis(:,ii),'Color', [0.870588 0.721569 0.529412],'Linewidth',0.1);
          hold on;
      end
  end
  hold off;
  
  drawnow;
  figxmin = min(min(XSurf));
  figxmax = max(max(XSurf));
  figymin = min(min(min(ZBED_xz)));
  figymax = max(max(max(ZBED_xz)));
  depthExaggerationTimes = 1;
  
  
  if NoAxisImages == 1
        figname = sprintf('%s\\y_%s.png', outputFolder_simulation, num2str(i));

        axis([figxmin figxmax figymin figymax]);
        
        xsize = (figxmax - figxmin)/depthExaggerationTimes; % using horizontal cross sections
        ysize = figymax - figymin;
        axis off;
        set(gca,'DataAspectRatio',[depthExaggerationTimes 1 1]);
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize/10 ysize/10]);
        set(gca, 'Position', get(gca, 'OuterPosition'));
        
        set(gcf,'GraphicsSmoothing','off');
        print('-dpng','-r100',figname);
   else
        
        axis([figxmin+5 figxmax-5 figymin figymax]);
        % axis off;
        box off;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);


        set(gca,'DataAspectRatio',[depthExaggerationTimes 1 1]);
%         set(gcf,'PaperUnits','inches','PaperPosition',[0 0 xsize ysize]);
        
        
        % xlabel('x','fontsize',10,'fontweight','normal','FontAngle','normal','fontname','Ariel');
        % ylabel('H','fontsize',10,'fontweight','normal','FontAngle','italic','fontname','Ariel');
        set(gca,'fontsize',10,'fontname','Ariel','XTickMode','manual','YTickMode','manual','LineWidth',0.5);
        
        figname = sprintf('%s\\A_%s.png', outputFolder_simulation, num2str(i));
        print('-dpng','-r600',figname);
        
        figname2 = sprintf('%s\\A_%s.eps', outputFolder_simulation, num2str(i));
        print(gcf,'-depsc','-vector',figname2); 
        epsclean(figname2);   
        
        
  end

  disp(num2str(i));

end