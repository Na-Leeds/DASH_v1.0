% DUNEINIT - A subroutine of NEWDUNES
% File written by Dave Rubin 
% Last modified in 3/2005

false = 0;
true = 1;
MORHRZ = true;
MORVRT = true;
FirstRun = true;

% GridSize = 100;
% SurfGridSpace = 0.5;
% EdgeGridSpace = 0.5;

gridnumber = GridSize/SurfGridSpace;
% N = -1;


%read input file or files
% eval(CurrentFig);
cd(casefolder);
% run('INPUTS_parameters.m');
run(parameterfilename);

CenterShift = round(GridSize/2);


[XSurf,YSurf] = meshgrid(-CenterShift+SurfGridSpace:SurfGridSpace:GridSize-CenterShift, -CenterShift+SurfGridSpace:SurfGridSpace:GridSize-CenterShift);
x = XSurf;
y = YSurf;
XEdge(1:GridSize/EdgeGridSpace) = -CenterShift+SurfGridSpace;
XEdge(GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1) = -CenterShift+SurfGridSpace:SurfGridSpace:GridSize-CenterShift;
XEdge(2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1) = GridSize-CenterShift;
XEdge(3*GridSize/EdgeGridSpace:4*GridSize/EdgeGridSpace-1) = GridSize-CenterShift:(-1)*SurfGridSpace:-CenterShift+SurfGridSpace;
YEdge(1:GridSize/EdgeGridSpace) = GridSize-CenterShift:(-1)*SurfGridSpace:-CenterShift+SurfGridSpace;
YEdge(GridSize/EdgeGridSpace:2*GridSize/EdgeGridSpace-1) = -CenterShift+SurfGridSpace;
YEdge(2*GridSize/EdgeGridSpace:3*GridSize/EdgeGridSpace-1) = -CenterShift+SurfGridSpace:SurfGridSpace:GridSize-CenterShift;
YEdge(3*GridSize/EdgeGridSpace:4*GridSize/EdgeGridSpace-1) = GridSize-CenterShift;

if TYPE == 3
  big = max(SPCNGF, SPCNGS);
  big = max(big, SPCNGT); % largest wavelenght in three sets
  BD = GridSize/big; % 16/04
end

% Prevent division by zero.
	if SPCNGF == 0; FD=10^5; end
	if SPCNGS == 0; SD=10^5; end
	if SPCNGT == 0; TD=10^5; end
	if SNSPF1 == 0; SNSPF1=10^5; end
	if SNSPF2 == 0; SNSPF2=10^5; end
	if SNSPS1 == 0; SNSPS1=10^5; end
	if SNSPS2 == 0; SNSPS2=10^5; end
	if SNSPT1 == 0; SNSPT1=10^5; end
	if SNSPT2 == 0 ; SNSPT2=10^5; end
	if SMPRDF == 0 ; SMPRDF=10^5; end
	if SMPRDS == 0 ; SMPRDS=10^5; end
	if SMPRDT == 0 ; SMPRDT=10^5; end
	if HTPRDF == 0 ; HTPRDF=10^5; end
	if HTPRDS == 0 ; HTPRDS=10^5; end
	if HTPRDT == 0 ; HTPRDT=10^5; end
	if VLPRDF == 0 ; VLPRDF=10^5; end
	if VLPRDS == 0 ; VLPRDS=10^5; end
	if VLPRDT == 0 ; VLPRDT=10^5; end
	if DEPPRD == 0 ; DEPPRD=10^5; end
	if SPCNGF ~= 0 ; FD=100/SPCNGF; end
	if SPCNGS ~= 0 ; SD=100/SPCNGS; end
	if SPCNGT ~= 0 ; TD=100/SPCNGT; end
    
% Make sure ZHORIZ is not equal to 1.0.
if abs(ZHORIZ-1.0) < 0.001 % Used to be 0.01
  ZHORIZ = 1.001;
end

	
% Calculate the migration direction of scour pits formed by intersecting
% bedform troughs of the first two sets of bedforms.  When specified in 
% the input paramaters, this calculation is used to rotate the bedforms 
% such that the sides of the block diagram are normal and parallel to the 
% axes of trough-shaped sets of cross-bedding..
if tan(TRENDF*2*pi/360) ~= tan(TRENDS*2*pi/360) % migration direction of bedform, first ~= second
  if (SPCNGF>0) && (SPCNGS>0)	% wavelenght of bedform > 0
    ANGLEA = atan2(VELOCF,tan((90+TRENDF-TRENDS)*2*pi/360)*VELOCF-VELOCS/sin((TRENDS-TRENDF) *2*pi /360)) *360/(2*pi)-90;
  end
end
if CHOICE ~= 0 
  ANGLEB = TRENDS-TRENDF;
  ANGLEC = TRENDT-TRENDF;
  TRENDF =- ANGLEA+(CHOICE-1)*90;
  TRENDS = TRENDF+ANGLEB;
  if SPCNGT ~= 0, TRENDT = TRENDF + ANGLEC; end
end
