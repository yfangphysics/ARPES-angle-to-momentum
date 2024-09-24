#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function YFATMmap()
	
	NewDataFolder /O root:Blue
	if (DataFolderExists("root:Blue:Map")==0)
		NewDataFolder root:Blue:Map
		YFMapInit()
	endif
	
	// Make the panel
	NewPanel /W=(300,100,1330,665) /K=1 /N=panelBlueMap as "DA30 Map Angle To Momentum"
	ModifyPanel cbRGB=(48896,52992,65280), fixedSize=1

   	// Make the spec list selection tools
	GroupBox specListGroupBox,pos={5,5},size={130,220},title="Cube Selection"
	NVAR selectedRow = root:Blue:Map:selectedRow
	ListBox specListBox,pos={10,25},size={120,170},frame=4,mode=1,proc=YFMapListBoxProc
	ListBox specListBox,listWave=root:Blue:Map:listWave,selRow=selectedRow
	ListBox specListBox help={"A text wave contain cube waves"}
	Button updateButton,pos={10,200},size={115,20},title="Update List",proc=YFMapButtonProc
	
	// Make the angle cube creation tools
	GroupBox angleGroupBox,pos={140,5},size={220,220},title="Angle Cube Creation"
	GroupBox angleGroupBox help={"Tools for displaying angle cube"}
	NVAR theta = root:Blue:Map:theta
	SetVariable thetaSetVar,pos={145,25},size={85,16},title="theta:",value=theta
	SetVariable thetaSetVar help={"The tilt value of the sample"}
	NVAR phid = root:Blue:Map:phid
	SetVariable phidSetVar,pos={235,25},size={85,16},title="phid:",value=phid
	SetVariable phidSetVar help={"The phi tilt value of the sample"}
	NVAR EF = root:Blue:Map:EF
	SetVariable EFSetVar,pos={145,50},size={170,16},title="Center Energy:",value=EF,proc=YFMapSetVarProc
	SetVariable EFSetVar help={"A special energy that will be in the center of one energy bin (it doesn't necessarily have to be the Fermi energy)"}
	NVAR EW = root:Blue:Map:EW
	SetVariable EWSetVar,pos={145,75},size={210,16},title="Energy Steps To Integrate:",value=EW,proc=YFMapSetVarProc
	SetVariable EWSetVar help={"The number of energy steps to include in each energy bin. This must be an odd number since it includes a central energy value and a few neighbors to either side."}
   NVAR rawDeltaE = root:Blue:Map:rawDeltaE
	TitleBox energy1TitleBox,pos={145,100},size={170,13},frame=0
	TitleBox energy2TitleBox,pos={145,125},size={170,13},frame=0
	Button angleButton,pos={145,150},size={165,20},title="Create Angle Cube",proc=YFMapButtonProc
	
	// Make the momentum cube creation tools
	GroupBox kGroupBox,pos={145,175},size={175,75},title="Momentum Cube Creation"
	GroupBox kGroupBox help={"Tools for converting an angle cube into a k cube whose axes are kx, ky, and energy"}
	NVAR pixels = root:Blue:Map:pixels
	SetVariable pixelsSetVar,pos={150,200},size={150,16},title="Pixels:",value=pixels,proc=YFMapSetVarProc,limits={2,inf,1}
	SetVariable pixelsSetVar help={"The largest dimension of the resulting k map, in pixels"}
	NVAR symmNeg = root:Blue:Map:symmNeg
	CheckBox checkNeg,pos={370,5},size={90,15},title="x/y --> -x/y",value= symmNeg,proc=YFMapCheckProc
	CheckBox checkNeg help={"Adds reflection symmetry about the x and y axes (four operations)"}
	NVAR symmXY = root:Blue:Map:symmXY
	CheckBox checkXY,pos={370,20},size={90,15},title="x --> y",value= symmXY,proc=YFMapCheckProc
	CheckBox checkXY help={"Adds reflection symmetry about the line y=x (two operations)"}
	NVAR symm3Rot = root:Blue:Map:symm3Rot
	CheckBox check3Rot,pos={460,5},size={70,15},title="3-fold rot.",value= symm3Rot,proc=YFMapCheckProc
	CheckBox check3Rot help={"Adds three-fold rotational symmetry (three operations)"}
	NVAR symm3Refl = root:Blue:Map:symm3Refl
	CheckBox check3Refl,pos={460,20},size={70,15},title="3-fold refl.",value= symm3Refl,proc=YFMapCheckProc
	CheckBox check3Refl help={"Adds three-fold rotational and reflection symmetry (six operations)"}
	NVAR rphi = root:Blue:Map:rphi
	SetVariable rphiSetVar,pos={370,35},size={85,16},title="rphi:",value=rphi
	SetVariable rphiSetVar help={"The rotation value of the map"}
	NVAR rotphi = root:Blue:Map:rotphi
	CheckBox checkrotphi,pos={460,35},size={70,15},title="rot. phi",value= rotphi,proc=YFMapCheckProc
	CheckBox checkrotphi help={"Adds phi rotation"}
	Button momentumButton,pos={370,60},size={150,20},title="Create Momentum Cube",proc=YFMapButtonProc
	
	// Make the miscellaneous tools
	GroupBox miscGroupBox,pos={365,80},size={175,145},title="Misc. Tools"
	Button stopKButton,pos={370,105},size={135,20},title="Stop Calculating K",proc=YFMapButtonProc
	Button stopKButton help={"Removes the k cube/map. Useful if you want to tweak the angle cube without automatically updating the k cube."}
	Button expAngleButton,pos={370,130},size={135,20},title="Export Angle Map",proc=YFMapButtonProc
	Button expAngleButton help={"Copies the angle map to the current folder"}
	Button expMomButton,pos={370,155},size={135,20},title="Export k Map",proc=YFMapButtonProc
	Button expMomButton help={"Copies the k-space map to the current folder"}
	Button expButton,pos={370,180},size={135,20},title="Export All",proc=YFMapButtonProc
	Button expButton help={"Exports the relevant cubes, maps, and settings to the current folder. It is smart enough to export the settings that were used to create the cubes, rather than the currently selected settings."}
	CheckBox checkNorm,pos={370,205},size={75,15},title="Normalize",value= 1,proc=YFMapCheckProc
	CheckBox checkNorm,help={"Decide if you want to normalize the data or not"}
	CheckBox checkNaN,pos={455,205},size={75,15},title="0 to NaN",value= 0,proc=YFMapCheckProc
	CheckBox checkNaN help={"Replaces all zeroes with NaNs in both maps (but not in the cubes themselves)"}

	// Make the energy slider
	GroupBox energyGroup,pos={5,230},size={540,35}
	GroupBox energyGroup help={"Use this slider to choose constant-energy slices of the angle and momentum cubes (referred to as angle maps and k maps)"}
	NVAR currentEnergy = root:Blue:Map:currentEnergy
	TitleBox energy3TitleBox,pos={10,240},size={120,12},frame=0
	TitleBox energy3TitleBox title="Energy: " + num2str(currentEnergy) + " eV"
	WAVE angleCube = root:Blue:Map:angleCube
	Slider energySlider,pos={125,240},size={415,16},value=currentEnergy,proc=YFMapSliderProc,vert= 0,ticks= 0
	Slider energySlider limits={DimOffset(angleCube,2),DimOffset(angleCube,2)+(DimSize(angleCube,2)-1)*DimDelta(angleCube,2),DimDelta(angleCube,2)}
	
	// Make the angle map display area
	GroupBox angleDispGroupBox,pos={5,270},size={540,270},title="Angle Map Display"
	Display/W=(10,290,540,535)/HOST=#
	RenameWindow #,GangleMap
	WAVE angleMap = root:Blue:Map:angleMap
	AppendImage angleMap
	ModifyImage angleMap ctab={*,*,ColdWarm,0}
	ModifyGraph wbRGB=(48896,52992,65280), gbRGB=(48896,52992,65280)
	ModifyGraph nticks=4,minor=1,fSize=12,standoff=0,tkLblRot(left)=0,btLen=3
	ModifyGraph margin(left)=50,margin(bottom)=35,margin(top)=10,margin(right)=25
	Label left "thetaa (deg)"
	Label bottom "phi (deg)"
	SetActiveSubwindow ##
	
	// Make the momentum map display area
	GroupBox momentumMapDispGroupBox,pos={550,5},size={470,535},title="Momentum Map Display"
	Display/W=(555,25,1005,535)/HOST=# 
	RenameWindow #,GkMap
	WAVE kMap = root:Blue:Map:kMap
	AppendImage kMap
	ModifyImage kMap ctab={*,*,ColdWarm,0}
	ModifyGraph wbRGB=(48896,52992,65280), gbRGB=(48896,52992,65280)
	//ModifyGraph nticks=4,minor=1,fSize=12,standoff=0,tkLblRot(left)=0,btLen=3
	//ModifyGraph margin(left)=50,margin(bottom)=35,margin(top)=5,margin(right)=5
	//ModifyGraph width=0, height={Aspect,DimDelta(kMap,1)*(DimSize(kMap,1)-1)/(DimDelta(kMap,0)*(DimSize(kMap,0)-1))}
	Label left "k\\By\\M (Å\\S-1\\M)"
	Label bottom "k\\Bx\\M (Å\\S-1\\M)"
	SetActiveSubwindow ##
	
	// Get ready to do stuff and update the display
	YFMapUpdateDisplay()
	
End

Function YFMapInit()
	
	// Store the current folder location
	String currDF = GetDataFolder(1)
	SetDataFolder root:Blue:Map
	
	// Declare global variables, strings, and waves
	Make /O /T /N=0 listWave
	Make /O /N=(2,2) angleMap
	Make /O /N=(2,2) kMap
	Make /O /N=(2,2,2) angleCube
	Make /O /N=(2,2,2) kCube
	Make /O /N=0 angleRecord
	Make /O /N=0 kRecord
	Variable /G selectedRow = -1
	Variable /G theta
	Variable /G phid
	Variable /G EF
	Variable /G EW 
	Variable /G status=0
	Variable /G doNaN = 0
	Variable /G doNorm = 1
	Variable /G currentEnergy
	Variable /G rawDeltaE
	Variable /G pixels = 200
	Variable /G symmNeg
	Variable /G symmXY
	Variable /G symm3Rot
	Variable /G symm3Refl
	Variable /G rotphi
	Variable /G rphi
	// Return to the current folder
	SetDataFolder currDF
	
	// Get the list of 2D text waves in the folder
	WAVE /T listWave = root:Blue:Map:listWave
	String list = WaveList("*",";","DIMS:3")
	Redimension /N=(ItemsInList(list)) listWave
	listWave = StringFromList(p,list)
	
End

Function YFMapUpdateDisplay()
	
	NVAR status = root:Blue:Map:status
	
	if (status == 0)	// No DA30 data has been chosen, so disable all the controls
		// Disable the angle cube controls
		SetVariable thetaSetVar disable=2
		SetVariable phidSetVar disable=2
		SetVariable EFSetVar disable=2
		SetVariable EWSetVar disable=2
		Button angleButton disable = 2
		// Disable the momentum cube controls
		SetVariable pixelsSetVar disable = 2
		SetVariable rphiSetVar disable=2
		Button momentumButton disable = 2
		CheckBox checkNeg disable = 2
		CheckBox checkXY disable = 2
		CheckBox check3Rot disable = 2
		CheckBox check3Refl disable = 2
		CheckBox checkrotphi disable = 2
		// Enable some of the misc. tools
		CheckBox checkNorm disable = 2
		CheckBox checkNaN disable = 2
		Button stopKButton disable = 2
		Button expAngleButton disable = 2
		Button expMomButton disable = 2
		Button expButton disable = 2
		// Disable the energy slider
		TitleBox energy3TitleBox title=""
		Slider energySlider disable=2
		// Hide the angle map display
		SetWindow panelBlueMap#GangleMap, hide=1
		// Hide the momentum map display
		SetWindow panelBlueMap#GkMap, hide=1
		
	elseif (status == 1)	// A DA30 data has been chosen, so enable the angle map controls
		// Enable the angle cube controls
		SetVariable thetaSetVar disable=0
		SetVariable phidSetVar disable=0
		SetVariable EFSetVar disable=0
		SetVariable EWSetVar disable=0
		Button angleButton disable = 0
		// Disable the momentum cube controls
		SetVariable pixelsSetVar disable = 2
		SetVariable rphiSetVar disable=2
		Button momentumButton disable = 2
		CheckBox checkNeg disable = 2
		CheckBox checkXY disable = 2
		CheckBox check3Rot disable = 2
		CheckBox check3Refl disable = 2
		CheckBox checkrotphi disable = 2
		// Enable some of the misc. tools
		CheckBox checkNorm disable = 0
		CheckBox checkNaN disable = 0
		Button stopKButton disable = 2
		Button expAngleButton disable = 2
		Button expMomButton disable = 2
		Button expButton disable = 2
		// Disable the energy slider
		TitleBox energy3TitleBox title=""
		Slider energySlider disable=2
		// Hide the angle map display
		SetWindow panelBlueMap#GangleMap, hide=1
		// Hide the momentum map display
		SetWindow panelBlueMap#GkMap, hide=1
		
	elseif (status == 2)	// An angle map has been created, so display it and enable the momentum map controls
		// Enable the angle cube controls
		SetVariable thetaSetVar disable=0
		SetVariable phidSetVar disable=0
		SetVariable EFSetVar disable=0
		SetVariable EWSetVar disable=0
		Button angleButton disable = 0
		// Enable the momentum cube controls
		SetVariable pixelsSetVar disable = 0
		SetVariable rphiSetVar disable=0
		CheckBox checkNeg disable = 0
		CheckBox checkXY disable = 0
		CheckBox check3Rot disable = 0
		CheckBox check3Refl disable = 0
		CheckBox checkrotphi disable = 0
		Button momentumButton disable = 0
		// Enable some of the misc. tools
		CheckBox checkNorm disable = 2
		CheckBox checkNaN disable = 0
		Button stopKButton disable = 2
		Button expAngleButton disable = 0
		Button expMomButton disable = 2
		Button expButton disable = 0
		// Enable the energy slider
		Slider energySlider disable=0
		// Show the angle map display
		SetWindow panelBlueMap#GangleMap, hide=0
		// Hide the momentum map display
		SetWindow panelBlueMap#GkMap, hide=1
		
	elseif (status == 3)	// A momentum map has been created, so display it and enable all the controls
		// Enable the angle cube controls
		SetVariable thetaSetVar disable=0
		SetVariable phidSetVar disable=0
		SetVariable EFSetVar disable=0
		SetVariable EWSetVar disable=0
		Button angleButton disable = 0
		// Enable the momentum cube controls
		SetVariable pixelsSetVar disable = 0
		SetVariable rphiSetVar disable=0
		CheckBox checkNeg disable = 0
		CheckBox checkXY disable = 0
		CheckBox check3Rot disable = 0
		CheckBox check3Refl disable = 0
		CheckBox checkrotphi disable = 0
		Button momentumButton disable = 0
		// Enable the misc. tools
		CheckBox checkNorm disable = 2
		CheckBox checkNaN disable = 2
		Button stopKButton disable = 0	
		Button expAngleButton disable = 0
		Button expMomButton disable = 0
		Button expButton disable = 0

		// Enable the energy slider
		Slider energySlider disable=0
		// Show the angle map display
		SetWindow panelBlueMap#GangleMap, hide=0
		// Hide the momentum map display
		SetWindow panelBlueMap#GkMap, hide=0
		
	endif
	
End

Function YFMapAngleCube()
   wave tempspec
	WAVE angleCube = root:Blue:Map:angleCube
	WAVE angleMap = root:Blue:Map:angleMap
	WAVE angleRecord = root:Blue:Map:angleRecord
	NVAR status = root:Blue:Map:status
	NVAR EF = root:Blue:Map:EF
	NVAR EW = root:Blue:Map:EW
	NVAR theta = root:Blue:Map:theta
	NVAR phid = root:Blue:Map:phid
	NVAR currentEnergy = root:Blue:Map:currentEnergy
	NVAR doNaN = root:Blue:Map:doNaN
	NVAR doNorm = root:Blue:Map:doNorm

	
	
	// Find the index of the closest energy value to EF
	
	Variable EFi = round((EF - DimOffset(tempspec,0))/DimDelta(tempspec,0))
	
	// Find the lowest index of any energy value we will use in the cube
	Variable Emini = EFi - 0.5*(EW-1) - EW*floor((EFi - 0.5*(EW-1))/EW)
	
	// Find the number of energy bins we will use in the cube
	Variable numE = floor((DimSize(tempSpec,0) - Emini)/EW)
	
	// Find the total number of q values
	Variable numQ = DimSize(tempSpec,1)
	
	// Prepare the cube
	redimension /N=(DimSize(tempspec,2), numQ, numE) angleCube
	SetScale /P x Dimoffset(tempspec,2), DimDelta(tempspec,2), angleCube
	SetScale /P y DimOffset(tempSpec,1), DimDelta(tempSpec,1), angleCube
	SetScale /P z DimOffset(tempSpec,0) + (Emini + 0.5*(EW-1))*DimDelta(tempSpec,0), EW*DimDelta(tempSpec,0), angleCube
	
	// Create the cube and the do normalization
	Variable i,j,k
	angleCube = 0
	for(j=0; j<numE; j+=1)
		for(k=0; k<EW; k+=1)
					angleCube[][][j] += tempSpec[Emini + EW*j + k][q][p]
		endfor
	endfor
	if(donorm)
	   variable rows,columns,tempcircle,tempaverage,conscircle
	   rows=DimSize(anglecube,0)
	   columns=DimSize(anglecube,1)
	   conscircle=0
	   do
	      tempcircle=0
	      do
		     make/n=(dimsize(anglecube,2))/o tewave=anglecube[tempcircle][conscircle][p]
		     Wavestats/Q tewave
		     tempaverage=V_avg
		     tewave=tewave/tempaverage
		     variable t
		     for(t=0;t<dimsize(anglecube,2);t++)
		         anglecube[tempcircle][conscircle][t]=tewave[t]
		     endfor
		     tempcircle+=1
		     KillWaves tewave								
	      while (tempcircle<rows)				
         conscircle+=1
       while (conscircle<columns)		
	endif
	
	// Update the angle map
	redimension /N=(DimSize(angleCube,0), DimSize(angleCube,1)) angleMap
	angleMap = angleCube[p][q][DimSize(angleCube,2)-1]
	SetScale /P x DimOffset(angleCube,0), DimDelta(angleCube,0), angleMap
	SetScale /P y DimOffset(angleCube,1), DimDelta(angleCube,1), angleMap
	
	// Remove zeros if desired
	if (doNan)
		angleMap =  (angleMap[p][q]==0) ? NaN : angleMap[p][q]
	endif
	
	// Update the energy slider
	Slider energySlider limits={DimOffset(angleCube,2),DimOffset(angleCube,2)+(DimSize(angleCube,2)-1)*DimDelta(angleCube,2),DimDelta(angleCube,2)}
	currentEnergy = DimOffset(angleCube,2)+(DimSize(angleCube,2)-1)*DimDelta(angleCube,2)
	Slider energySlider value=currentEnergy
	TitleBox energy3TitleBox title="Energy: " + num2str(currentEnergy) + " eV"
	
	// Update the status
	if (status < 2)
		status = 2
	endif
	
	// Store a record of the settings
	Redimension /N=6 angleRecord
	angleRecord[0] = theta
	angleRecord[1] = phid
	angleRecord[2] = EF
	angleRecord[3] = EW
	angleRecord[4] = doNorm
	angleRecord[5] = doNaN
	
	killwaves tempspec
End



Function YFMapKCube()
   	WAVE angleCube = root:Blue:Map:angleCube
	WAVE kCube = root:Blue:Map:kCube
	WAVE kMap = root:Blue:Map:kMap
	WAVE kRecord = root:Blue:Map:kRecord
	NVAR pixels = root:Blue:Map:pixels
	NVAR currentEnergy = root:Blue:Map:currentEnergy
	NVAR doNaN = root:Blue:Map:doNaN
	NVAR status = root:Blue:Map:status
	NVAR theta = root:Blue:Map:theta
	NVAR phid = root:Blue:Map:phid
	
	
	// Confirm that we should be doing this
	if (status < 2)
		return 0
	endif

   Make/d /O /N=(DimSize(angleCube,0), DimSize(angleCube,1)) angleMap
	SetScale /P x DimOffset(angleCube,0), DimDelta(angleCube,0), angleMap
	SetScale /P y DimOffset(angleCube,1), DimDelta(angleCube,1), angleMap
	angleMap = angleCube[p][q][DimSize(angleCube,2)-1]
	
   MakeList(angleMap, DimOffset(angleCube,2) + (DimSize(angleCube,2)-1)*DimDelta(angleCube,2))
   wave momentumlist
	// Now find the k-space bounding box of the map
	Variable kxMin, kxMax, kyMin, kyMax
	Make /O /N=(DimSize(momentumList,0)) tempWave=0
	tempWave = momentumList[p][0]
	WaveStats /Q tempWave
	kxMin = V_min
	kxMax = V_max
	tempWave = momentumList[p][1]
	WaveStats /Q tempWave
	kyMin = V_min
	kyMax = V_max
	
	// Determine the final size for the map
	Variable ratio = (kxMax - kxMin)/(kyMax - kyMin)
	Variable mapSizeX, mapSizeY, delta
	if (ratio > 1)
		mapSizeX = pixels
		mapSizeY = floor(pixels/ratio)
		delta = (kxMax - kxMin)/(pixels - 1)
	else
		mapSizeX = floor(pixels*ratio)
		mapSizeY = pixels
		delta = (kyMax - kyMin)/(pixels - 1)
	endif
	
	// Set up the k cube and the Jacobian
	Redimension /N=(mapSizeX, mapSizeY, DimSize(anglecube,2)) kCube
	SetScale /P x kxMin, delta, kCube
	SetScale /P y kyMin, delta, kCube
	SetScale /P z DimOffset(anglecube,2), DimDelta(anglecube,2), kCube
	kCube = 0
	Make /O /N=(mapSizeX, mapSizeY) BlueJacobian
	BlueJacobian = 0
	
	// Make the k cube
	Variable pee, qew
	Variable i,j,u,v
	Make /O /N=(DimSize(angleCube,0)*ceil(DimDelta(angleCube,0)/DimDelta(angleCube,1)), DimSize(angleCube,1)) YFUpsampled
	SetScale /I x DimOffset(angleCube,0), DimOffset(angleCube,0) + (DimSize(angleCube,0)-1)*DimDelta(angleCube,0), YFUpsampled
	SetScale /P y DimOffset(angleCube,1), DimDelta(angleCube,1), YFUpsampled
   for(i=0;i<DimSize(anglecube,2);i++)
       angleMap = angleCube[p][q][i]
       // Up-scale the angle map to compensate for large phi step sizes
		YFUpsampled = Interp2D(angleMap,x,y)
       
       MakeList(YFUpsampled, DimOffset(angleCube,2) + i*DimDelta(angleCube,2))
		// Store the momentum list in the k cube
		BlueJacobian = 0
		for (j=0; j<DimSize(momentumList,0); j+=1)
			pee = round((momentumList[j][0] - DimOffset(kCube,0))/DimDelta(kCube,0))
			qew = round((momentumList[j][1] - DimOffset(kCube,1))/DimDelta(kCube,1))
			if (pee >= 0 && pee < DimSize(kCube,0) && qew >= 0 && qew < DimSize(kCube,1))
				kCube[pee][qew][i] += momentumList[j][2]
				BlueJacobian[pee][qew] += 1
			endif
		endfor
		for(u=0;u<DimSize(kCube,0);u+=1)
		    for(v=0;v<DimSize(kCube,1);v+=1)
		        kCube[u][v][i] /= max(1,BlueJacobian[u][v])
		    endfor
		endfor
	endfor
	
	// Clean up
	KillWaves momentumList, BlueJacobian, tempWave,angleMap, YFUpsampled
	
	// Update the k map and the graph
	Redimension /N=(DimSize(kCube,0), DimSize(kCube,1)) kMap
	kMap = kCube[p][q][round((currentEnergy - DimOffset(kCube,2))/DimDelta(kCube,2))]
	SetScale /P x DimOffset(kCube,0), DimDelta(kCube,0), kMap
	SetScale /P y DimOffset(kCube,1), DimDelta(kCube,1), kMap
	ModifyGraph /W=PanelBlueMap#GkMap width=0, height={Aspect,1/ratio}
	
	// Remove zeros if desired
	if (doNan)
		kMap =  (kMap[p][q]==0) ? NaN : kMap[p][q]
	endif
	
	// Update the status if necessary
	if (status < 3)
		status = 3
	endif
	
	// Make a record of the settings used to create the cube
	Redimension /N=1 kRecord
	kRecord[0] = donan
End


function MakeList(angleMap, energy)
  	
	WAVE angleMap
	Variable energy
	NVAR theta = root:Blue:Map:theta
	NVAR phid = root:Blue:Map:phid
	NVAR symmNeg = root:Blue:Map:symmNeg
	NVAR symmXY = root:Blue:Map:symmXY
	NVAR symm3Rot = root:Blue:Map:symm3Rot
	NVAR symm3Refl = root:Blue:Map:symm3Refl
	NVAR rotphi=root:blue:map:rotphi
	NVAR rphi = root:Blue:Map:rphi
	
	Make /O /N=(DimSize(angleMap,0)*DimSize(angleMap,1),3) momentumList

   variable i,j,n=0
   variable kx,ky
   variable phi,thetaa
	// Make a list and convert to momentum
	for (i=0; i<DimSize(angleMap,0); i+=1)
	   phi=Dimoffset(angleMap,0)+i*Dimdelta(angleMap,0)
		for (j=0; j<DimSize(angleMap,1); j+=1)
		   
		   thetaa=Dimoffset(angleMap,1)+j*Dimdelta(angleMap,1)
		   
		   kx=0.512*sqrt(energy)*sin((phi+phid)*pi/180)*cos(thetaa*pi/180)
		   ky=0.512*sqrt(energy)*(cos(theta*pi/180)*sin(thetaa*pi/180)-sin(theta*pi/180)*cos((phi+phid)*pi/180)*cos(thetaa*pi/180))

			momentumList[n][0] = kx
			momentumList[n][1] = ky
			momentumList[n][2] = angleMap[i][j]
			n += 1
		endfor
	endfor
	
	// Apply symmetry
	Variable num
	if (symmNeg)
		num = DimSize(momentumList,0)
		Redimension /N=(num*4,3) momentumList
		
		momentumList[num, 2*num-1][0] = -momentumList[p-n][0]
		momentumList[num, 2*num-1][1] = momentumList[p-n][1]
		momentumList[num, 2*num-1][2] = momentumList[p-n][2]
		
		momentumList[2*num, 3*num-1][0] = momentumList[p-2*num][0]
		momentumList[2*num, 3*num-1][1] = -momentumList[p-2*num][1]
		momentumList[2*num, 3*num-1][2] = momentumList[p-2*num][2]
		
		momentumList[3*num, 4*num-1][0] = -momentumList[p-3*num][0]
		momentumList[3*num, 4*num-1][1] = -momentumList[p-3*num][1]
		momentumList[3*num, 4*num-1][2] = momentumList[p-3*num][2]
	endif
	
	if (symmXY)
		num = DimSize(momentumList,0)
		Redimension /N=(num*2,3) momentumList
		
		momentumList[num, 2*num-1][0] = momentumList[p-num][1]
		momentumList[num, 2*num-1][1] = momentumList[p-num][0]
		momentumList[num, 2*num-1][2] = momentumList[p-num][2]
	endif
	
	
	if (symm3Rot)
		num = DimSize(momentumList,0)
		Redimension /N=(num*3,3) momentumList
		
		momentumList[num,2*num-1][0] = -(1/2)*momentumList[p-num][0] - (sqrt(3)/2)*momentumList[p-num][1]
		momentumList[num,2*num-1][1] = (sqrt(3)/2)*momentumList[p-num][0] - (1/2)*momentumList[p-num][1]
		momentumList[num,2*num-1][2] = momentumList[p-num][2]
		
		momentumList[2*num,3*num-1][0] = -(1/2)*momentumList[p-2*num][0] + (sqrt(3)/2)*momentumList[p-2*num][1]
		momentumList[2*num,3*num-1][1] = -(sqrt(3)/2)*momentumList[p-2*num][0] - (1/2)*momentumList[p-2*num][1]
		momentumList[2*num,3*num-1][2] = momentumList[p-2*num][2]
	endif
	
	if (rotphi)
	   num = DimSize(momentumList,0)
		Redimension /N=(num*3,3) momentumList
		
		momentumList[num,2*num-1][0] = -(1/2)*momentumList[p-num][0] - (sqrt(3)/2)*momentumList[p-num][1]
		momentumList[num,2*num-1][1] = (sqrt(3)/2)*momentumList[p-num][0] - (1/2)*momentumList[p-num][1]
		momentumList[num,2*num-1][2] = momentumList[p-num][2]
		
		momentumList[2*num,3*num-1][0] = -(1/2)*momentumList[p-2*num][0] + (sqrt(3)/2)*momentumList[p-2*num][1]
		momentumList[2*num,3*num-1][1] = -(sqrt(3)/2)*momentumList[p-2*num][0] - (1/2)*momentumList[p-2*num][1]
		momentumList[2*num,3*num-1][2] = momentumList[p-2*num][2]
		
	   num = DimSize(momentumList,0)
	   duplicate/o momentumlist,ttc
	   variable ccc
	   for(ccc=0;ccc<num;ccc++)
	       momentumlist[ccc][0]=cos(rphi/180*pi)*ttc[ccc][0]-sin(rphi/180*pi)*ttc[ccc][1]
	       momentumlist[ccc][1]=sin(rphi/180*pi)*ttc[ccc][0]+cos(rphi/180*pi)*ttc[ccc][1]
	       momentumlist[ccc][2]=ttc[ccc][2]
	   endfor
	   killwaves ttc
	endif
	   
	
	if (symm3Refl)
		num = DimSize(momentumList,0)
		Redimension /N=(num*6,3) momentumList
		
		// Rotate by 120 degrees
		momentumList[num,2*num-1][0] = -(1/2)*momentumList[p-num][0] - (sqrt(3)/2)*momentumList[p-num][1]
		momentumList[num,2*num-1][1] = (sqrt(3)/2)*momentumList[p-num][0] - (1/2)*momentumList[p-num][1]
		momentumList[num,2*num-1][2] = momentumList[p-num][2]
		
		// Rotate by 240 degrees
		momentumList[2*num,3*num-1][0] = -(1/2)*momentumList[p-2*num][0] + (sqrt(3)/2)*momentumList[p-2*num][1]
		momentumList[2*num,3*num-1][1] = -(sqrt(3)/2)*momentumList[p-2*num][0] - (1/2)*momentumList[p-2*num][1]
		momentumList[2*num,3*num-1][2] = momentumList[p-2*num][2]
		
		// Flip about the y axis
		momentumList[3*num,4*num-1][0] = -momentumList[p-3*num][0]
		momentumList[3*num,4*num-1][1] = momentumList[p-3*num][1]
		momentumList[3*num,4*num-1][2] = momentumList[p-3*num][2]
		
		// Flip about the y axis, rotate by 120 degrees
		momentumList[4*num,5*num-1][0] = (1/2)*momentumList[p-4*num][0] - (sqrt(3)/2)*momentumList[p-4*num][1]
		momentumList[4*num,5*num-1][1] = -(sqrt(3)/2)*momentumList[p-4*num][0] - (1/2)*momentumList[p-4*num][1]
		momentumList[4*num,5*num-1][2] = momentumList[p-4*num][2]
		
		// Flip about the y axis, rotate by 240 degrees
		momentumList[5*num,6*num-1][0] = (1/2)*momentumList[p-5*num][0] + (sqrt(3)/2)*momentumList[p-5*num][1]
		momentumList[5*num,6*num-1][1] = (sqrt(3)/2)*momentumList[p-5*num][0] - (1/2)*momentumList[p-5*num][1]
		momentumList[5*num,6*num-1][2] = momentumList[p-5*num][2]
	endif
end

Function YFMapSetVarProc(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	
	NVAR EW = root:Blue:Map:EW
	NVAR EF = root:Blue:Map:EF
	NVAR rawDeltaE = root:Blue:Map:rawDeltaE
	NVAR pixels = root:Blue:Map:pixels
   WAVE tempspec
   
	strswitch(ctrlName)
		case "EWSetVar":
			// Convert to an allowed value
			EW = max(1,1 + 2*round((varNum-1)/2))
			
			// Check to see if the value is too big
			Variable EFi = DimOffset(tempSpec,0) + DimDelta(tempSpec,0)*round((EF - DimOffset(tempSpec,0))/DimDelta(tempSpec,0))
			Variable diff = min(EFi - DimOffset(tempSpec,0), DimOffset(tempSpec,0) + (DimSize(tempSpec,0)-1)*DimDelta(tempSpec,0) - EFi)
			if (0.5*(EW-1)*rawDeltaE > diff)
				// The value is too big, so revert back to 1
				EW = 1
			endif
			
			// Update the text display
			TitleBox energy2TitleBox title="Your Chosen Bin Width: " + num2str(rawDeltaE*(EW-1)) + " eV"
			break
		case "pixelsSetVar":
			// Make sure the value is an integer
			pixels = round(varNum)
			break

	endswitch
End

Function YFMapSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	strswitch(sa.ctrlName)
		case "energySlider":
			if(sa.eventCode==3||sa.eventCode==9)
				NVAR currentEnergy = root:Blue:Map:currentEnergy
				NVAR status = root:Blue:Map:status
			       NVAR doNaN = root:Blue:Map:doNaN
				WAVE angleCube = root:Blue:Map:angleCube
				WAVE angleMap = root:Blue:Map:angleMap
				WAVE kCube = root:Blue:Map:kCube
				WAVE kMap = root:Blue:Map:kMap
				
				// Update the current energy
				currentEnergy = sa.curval
				
				// Update the energy display
				TitleBox energy3TitleBox title="Energy: " + num2str(currentEnergy) + " eV"
				
				
				// Update the angleMap
				angleMap = angleCube[p][q][round((currentEnergy-DimOffset(angleCube,2))/DimDelta(angleCube,2))]
				
				Redimension /N=(DimSize(kCube,0), DimSize(kCube,1)) kMap
				SetScale /P x DimOffset(kCube,0), DimDelta(kCube,0), kMap
	          SetScale /P y DimOffset(kCube,1), DimDelta(kCube,1), kMap
	          if (status == 3)
				    // Update the kMap
				    kMap = kCube[p][q][round((currentEnergy-DimOffset(kCube,2))/DimDelta(kCube,2))]
				    	// Remove zeros if desired
	                       if (doNan)
		                      kMap =  (kMap[p][q]==0) ? NaN : kMap[p][q]
	                       endif
				endif
				
				
			endif			
			break
	endswitch
	
End

Function YFMapButtonProc(ctrlName)
	String ctrlName
	
	WAVE /T listWave = root:Blue:Map:listWave
	NVAR stage = root:Blue:Map:status
	NVAR doNorm = root:Blue:Map:doNorm
	NVAR doNan = root:Blue:Map:doNan
	NVAR selectedRow = root:Blue:Map:selectedRow
	
	strswitch(ctrlName)
		case "updateButton":
			// Update the wave list
			string list = WaveList("*",";","DIMS:3")
			Redimension /N=(ItemsInList(list)) listWave
			listWave = StringFromList(p,list)
			
			// Reset the status and update the display
			selectedRow = -1
			ListBox specListBox,selRow = selectedRow
			stage = 0
			donorm=1
			CheckBox checkNorm value = 1
			donan=0
			CheckBox checkNan value = 0
			YFMapUpdateDisplay()
			break
		case "angleButton":
			YFMapAngleCube()
			if (stage == 3)
				YFMapKCube()
			endif
			YFMapUpdateDisplay()
			break

		case "momentumButton":
			YFMapKCube()
			YFMapUpdateDisplay()
			break
		case "expMomButton":
			Duplicate /O root:Blue:Map:kMap kMap
			break
		case "expAngleButton":
			Duplicate /O root:Blue:Map:angleMap angleMap
			break
		case "stopKButton":
			stage = 2
			YFMapUpdateDisplay()
			break
		case "expButton":
			WAVE angleRecord = root:Blue:Map:angleRecord
			WAVE kRecord = root:Blue:Map:kRecord
			WAVE angleCube0 = root:Blue:Map:angleCube
			WAVE kCube0 = root:Blue:Map:kCube
			WAVE angleMap0 = root:Blue:Map:angleMap
			WAVE kMap0 = root:Blue:Map:kMap
			NVAR currentEnergy0 = root:Blue:Map:currentEnergy
			
			if (stage <= 1)	// This should never happen, but if it does, don't do anything
				return 0
			endif
			
			if (stage > 1)	// There is an angle cube
				Duplicate /O angleCube0 angleCube
				Duplicate /O angleMap0 angleMap
				Variable /G theta = angleRecord[0]
				Variable /G phid = angleRecord[1]
				Variable /G EF = angleRecord[2]
				Variable /G EW = angleRecord[3]
				Variable /G doNorm = angleRecord[4]
				Variable /G doNan = angleRecord[5]
			endif
			
			if (stage == 3)	// There is also a k map
				Duplicate /O kCube0 kCube
				Duplicate /O kMap0 kMap
				doNan = kRecord[0]
			endif
			Variable /G status = stage
			Variable /G currentEnergy = currentEnergy0
			break
	endswitch
End

Function YFMapListBoxProc(ctrlName,row,col,event) : ListBoxControl
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	NVAR status = root:Blue:Map:status
	NVAR EF = root:Blue:Map:EF
	NVAR rawDeltaE = root:Blue:Map:rawDeltaE
	NVAR EW = root:Blue:Map:EW
	WAVE /T specList = root:Blue:Map:specList
	
	if (event == 3)
		WAVE /T listWave = root:Blue:Map:listWave
		NVAR selectedRow = root:Blue:Map:selectedRow
		if (selectedRow != row)
			selectedRow = row
			if (selectedRow > numpnts(listWave) - 1)
				selectedRow = -1
				ListBox specListBox,selRow=selectedRow
			else
			   if(status==0)
	             status=1
	          endif
	       endif			
			YFMapUpdateDisplay()
		endif
		duplicate/o $listwave[selectedRow],tempspec
	   EF = DimOffset(tempSpec,0) + (DimSize(tempSpec,0)-1)*DimDelta(tempSpec,0)
		rawDeltaE = DimDelta(tempSpec,0)
		// Update some displays
		TitleBox energy1TitleBox title="Spectrum Energy Step: " + num2str(rawDeltaE) + " eV"
		TitleBox energy2TitleBox title="Your Chosen Bin Width: " + num2str(rawDeltaE*(EW-1)) + " eV"
	endif
	return 0
End

Function YFMapCheckProc(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked

	
	strswitch (ctrlName)
		case "checkNorm":
			NVAR doNorm = root:Blue:Map:doNorm
			if (checked)
		      donorm=1
             CheckBox checkNorm value = 1
			else
			   donorm=0
             CheckBox checkNorm value = 0
			endif
			break
		case "checkNaN":
			NVAR doNaN = root:Blue:Map:doNaN
			if (checked)
			   donan=1
             CheckBox checkNan value = 1				
			else
			   donan=0
             CheckBox checkNan value = 0				
			endif
			break
		case "checkNeg":
			NVAR symmNeg = root:Blue:Map:symmNeg
			if (checked)
				symmNeg = 1
			else
				symmNeg = 0
			endif
			break
		case "checkXY":
			NVAR symmXY = root:Blue:Map:symmXY
			if (checked)
				symmXY = 1
			else
				symmXY = 0
			endif
			break
		case "check3Rot":
			NVAR symm3Rot = root:Blue:Map:symm3Rot
			NVAR symm3Refl = root:Blue:Map:symm3Refl
			if (checked)
				symm3Rot = 1
				symm3Refl = 0
				CheckBox check3Refl value = 0
			else
				symm3Rot = 0
			endif
			break
		case "check3Refl":
			NVAR symm3Refl = root:Blue:Map:symm3Refl
			NVAR symm3Rot = root:Blue:Map:symm3Rot
			if (checked)
				symm3Refl = 1
				symm3Rot = 0
				CheckBox check3Rot value = 0
			else
				symm3Refl = 0
			endif
			break
	  case "checkrotphi":
			NVAR rotphi = root:Blue:Map:rotphi
			if (checked)
				rotphi = 1
			else
				rotphi = 0
			endif
			break
	endswitch
	
End
