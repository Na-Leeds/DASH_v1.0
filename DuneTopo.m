% DUNETOPO - A subroutine of NEWDUNES that calculates the topography of bedforms at each time step
% File written by Dave Rubin 
% Last modified in 3/2005

TIME_temp = TIME_temp + 1;
TIME_list(TIME_temp,1) = TIME;


% Calculate deposition.
DPOSIT = (DEPRAT*TIME)+DEPCHG*(-DEPPRD/(2*pi))*cos(((DEPFAZ-90)*2*pi/360)+(TIME*2*pi/DEPPRD));
DPOSIT_temp = DPOSIT_temp + 1;
DPOSIT_List(DPOSIT_temp,1) = DPOSIT;
% DEPRAT: Rate of deposition, lenght/time
% DEPCHG: Amplitude of cycle in rate of deposition, length/time
% DEPPRD: Period of cycle in rate of deposition, arbitrary units of time
% DEPFAZ: Phase of cycle in rate of deposition (in degrees); if 0, the rate of deposition is at a minimum at time 0





% Calculate bedform asymmetry.
PROFF = ((1-SMTRYF)+SMCHGF*sin(((-TIME/SMPRDF)-SMFAZF/360)*pi*2.0))*pi/2.0;
PROFS = ((1-SMTRYS)+SMCHGS*sin(((-TIME/SMPRDS)-SMFAZS/360)*pi*2.0))*pi/2.0;
PROFT = ((1-SMTRYT)+SMCHGT*sin(((-TIME/SMPRDT)-SMFAZT/360)*pi*2.0))*pi/2.0;

% Calculate bedform sizes; steepness
sizef = HTRTOF+(HTCHGF*sin(((TIME/HTPRDF)*2*pi)+((HTFAZF-90)*2*pi/360.0)));
sizes = HTRTOS+(HTCHGS*sin(((TIME/HTPRDS)*2*pi)+((HTFAZS-90)*2*pi/360.0)));
sizet = HTRTOT+(HTCHGT*sin(((TIME/HTPRDT)*2*pi)+((HTFAZT-90)*2*pi/360.0)));

% sizeIncreaseRate = -0.005; % positive increase in size; negative decrease in size
sizef = sizef - sizef * sizeIncreaseRate * (TIMEEnd - TIME + 1);
sizes = sizes - sizes * sizeIncreaseRate * (TIMEEnd - TIME + 1);
sizet = sizet - sizet * sizeIncreaseRate * (TIMEEnd - TIME + 1);

sizef_temp = sizef_temp + 1;
sizef_List(sizef_temp,1) = sizef;

% Calculate bedform locations,displacement of the sine curves simulates bedform migration
dispfd = VELOCF*TIME+VLCHGF*(-VLPRDF/(2*pi))*cos(((VLFAZF-90)*2*pi/360)+(TIME*2*pi/VLPRDF));
dispsd = VELOCS*TIME+VLCHGS*(-VLPRDS/(2*pi))*cos(((VLFAZS-90)*2*pi/360)+(TIME*2*pi/VLPRDS));
disptd = VELOCT*TIME+VLCHGT*(-VLPRDT/(2*pi))*cos(((VLFAZT-90)*2*pi/360)+(TIME*2*pi/VLPRDT));
% VELOCF: (22) Mean migration speed of bedform
% VLCHGF: (23) Amplitude of speed cycle
% VLPRDF: (24) Period of speed cycle
% VLFAZF: (25) Phase of speed cycle,in degrees

dispfd_temp = dispfd_temp + 1;
dispfd_List(dispfd_temp,1) = dispfd;



yfd = x*sin(TRENDF*2*pi/360)+y*cos(TRENDF*2*pi/360);
ysd = x*sin(TRENDS*2*pi/360)+y*cos(TRENDS*2*pi/360);
ytd = x*sin(TRENDT*2*pi/360)+y*cos(TRENDT*2*pi/360);

xfd = (x*cos(TRENDF*2*pi/360)-y*sin(TRENDF*2*pi/360)-dispfd-SPCNGF*PHASEF/360)-...
(SNMGF1*sin((yfd*2*pi/SNSPF1)+((SNFZF1*2*pi/360)+(TIME*SNVLF1*2*pi/SNSPF1))))-...
(SNMGF2*sin((yfd*2*pi/SNSPF2)+((SNFZF2*2*pi/360)+(TIME*SNVLF2*2*pi/SNSPF2))));

xsd = (x*cos(TRENDS*2*pi/360)-y*sin(TRENDS*2*pi/360)-...
      dispsd-SPCNGS*PHASES/360)-(SNMGS1*sin((ysd*2*pi/SNSPS1)+...
      ((SNFZS1*2*pi/360)+(TIME*SNVLS1*2*pi/SNSPS1))))-...
      (SNMGS2*sin((ysd*2*pi/SNSPS2)+((SNFZS2*2*pi/360)+...
      (TIME*SNVLS2*2*pi/SNSPS2))));
xtd = (x*cos(TRENDT*2*pi/360)-y*sin(TRENDT*2*pi/360)-...
      disptd-SPCNGT*PHASET/360)-(SNMGT1*sin((ytd*2*pi/SNSPT1)+...
      ((SNFZT1*2*pi/360)+(TIME*SNVLT1*2*pi/SNSPT1))))-...
      (SNMGT2*sin((ytd*2*pi/SNSPT2)+((SNFZT2*2*pi/360)+...
      (TIME*SNVLT2*2*pi/SNSPT2))));
  
% check
% figure; contour(xfd); figure; contour(yfd);

 
% Superimpose sets of dunes
if TYPE == 1 % simple addition
  zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD);
  shape = 1;
  zfd = zfd*sizef;
  zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD)*shape*sizes;
  ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD)*shape*sizet;
  z = zfd + zsd + ztd + DPOSIT;
  z = max (z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*ELVMIN+DPOSIT); % ELVMIN: elevation of interdune flats
  
elseif TYPE == 2  % proportional to the main bedform
  zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD);
  shape = ((7.5/FD)-zfd)/(15/FD) ;
  zfd = zfd*sizef;
  zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD).*shape*sizes;
  ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD).*shape*sizet;
  z = zfd + zsd + ztd + DPOSIT;
  z = max (z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*ELVMIN+DPOSIT);

elseif TYPE == 3 % locally highest
  zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD+7.5/FD)*sizef;
  zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD+7.5/SD)*sizes;
  ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD+7.5/TD)*sizet;
  z = max(zfd, zsd); z = max(z, ztd);
  z = z + DPOSIT;
  z = max (z, (7.5*(1+ELVMIN)/BD+DPOSIT) );
  
elseif TYPE == 4 % inversely proportional
  zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD);
  shape = 1-((7.5/FD)-zfd)/(15/FD);
  zfd = zfd*sizef;
  zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD).*shape*sizes;
  ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD).*shape*sizet;
  z = zfd + zsd + ztd + DPOSIT;
  z = max (z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*ELVMIN+DPOSIT);
  
elseif TYPE == 5
  zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD)*sizef;
  zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD+7.5/SD)*sizes;
  ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD+7.5/TD)*sizet;
  z = zfd + max(zsd, ztd) + DPOSIT;
  z = max (z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*ELVMIN+DPOSIT);
  
elseif TYPE == 6 % there is error here
  zfd = (-6*sin(xfd*FD*pi/50)/FD-1.5*sin((xfd*FD*pi/25)+PROFF)/FD);
  shape = ((7.5/FD)-zfd)/(15/FD);
  zfd = zfd*sizef;
  zsd = (-6*sin(xsd*SD*pi/50)/SD-1.5*sin((xsd*SD*pi/25)+PROFS)/SD).*shape*sizes;
  ztd = (-6*sin(xtd*TD*pi/50)/TD-1.5*sin((xtd*TD*pi/25)+PROFT)/TD).*shape*sizet;
  z = zfd + zsd + ztd + DPOSIT;
  z = max (z, ((7.5/FD)+(7.5/SD)+(7.5/TD))*ELVMIN+DPOSIT);
end

% Determine high and low points on bedform surface, and set initial values
% of elevation arrays (zcont and ZBED).
% ZBED: z of the domain edge, from bottomleft clockwise direction
if FirstRun
  zmin = min(min(z));
  zmax = max(max(z));
  zref = zmin + (zmax-zmin)*ZHORIZ;
  if(zref < -30.0); zref = -30.0; end
  % ZBED(GridSize/EdgeGridSpace:-EdgeGridSpace:1) = z(:,1);
  ZBED(GridSize/EdgeGridSpace:-1:1) = z(:,1);
  ZBED(GridSize/EdgeGridSpace:-1:1) = z(:,1);
  ZBED(GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1) = z(1,:);
  ZBED(2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1) = z(:,100);  %%
  % ZBED(4*GridSize/EdgeGridSpace-1:-EdgeGridSpace:3*GridSize/EdgeGridSpace) = z(100,:);  %%
  ZBED(4*GridSize/EdgeGridSpace-1:-1:3*GridSize/EdgeGridSpace) = z(100,:); 
  ZBED = min(ZBED,zref);
  ZBED = max(ZBED,-30);
  
  
  
end

%Define ZBED
if size(z,1) > 1
  %zcont = min(z,zcont);
  ZBED(GridSize/EdgeGridSpace:-1:1) = min(ZBED(GridSize/EdgeGridSpace:-1:1),z(:,1)');
  ZBED(GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1) = min(ZBED(GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1),z(1,:));
  ZBED(2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1) = min(ZBED(2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1),z(:,gridnumber)');
  ZBED(4*GridSize/EdgeGridSpace-1:-1:3*GridSize/EdgeGridSpace) = min(ZBED(4*GridSize/EdgeGridSpace-1:-1:3*GridSize/EdgeGridSpace),z(gridnumber,:));   %%
  %  ZBED = min(ZBED,zref); 
  ZBED = max(ZBED,-30);
else
    

  ZBED = min(z, ZBED); % ZBED is pervious line
  ZBED = max(ZBED,-30);
  
end



