#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function YFATMSpec()
	
	NewDataFolder /O root:Blue
	if (DataFolderExists("root:Blue:Map")==0)
		NewDataFolder root:Blue:Map
		YFSpecInit()
	endif
	
	// Make the panel
	NewPanel /W=(300,100,1130,665) /K=1 /N=panelBlueMap as "Cut Angle To Momentum"
	ModifyPanel cbRGB=(48896,52992,65280), fixedSize=1

   	// Make the spec list selection tools
	GroupBox specListGroupBox,pos={5,5},size={130,220},title="Spec Selection"
	NVAR selectedRow = root:Blue:Map:selectedRow
	ListBox specListBox,pos={10,25},size={120,170},frame=4,mode=1,proc=YFSpecListBoxProc
	ListBox specListBox,listWave=root:Blue:Map:listWave,selRow=selectedRow
	ListBox specListBox help={"A text wave contain spec waves"}
	Button updateButton,pos={10,200},size={115,20},title="Update List",proc=YFSpecButtonProc
	
	// Make the momentum cut creation tools
	NVAR theta = root:Blue:Map:theta
	SetVariable thetaSetVar,pos={145,25},size={85,16},title="theta:",value=theta
	SetVariable thetaSetVar help={"The tilt value of the sample"}
	NVAR phid = root:Blue:Map:phid
	SetVariable phidSetVar,pos={235,25},size={85,16},title="phid:",value=phid
	SetVariable phidSetVar help={"The phi tilt value of the sample"}
	NVAR EF = root:Blue:Map:EF
	SetVariable EFSetVar,pos={145,50},size={85,16},title="EF:",value=EF
	SetVariable EFSetVar help={"A special energy that will be in the center of one energy bin (it doesn't necessarily have to be the Fermi energy)"}
	NVAR bin = root:Blue:Map:bin
	SetVariable binSetVar,pos={235,50},size={85,16},title="bin:",value=bin
	SetVariable binSetVar help={"bin"}
	
	// Make the miscellaneous tools
	GroupBox miscGroupBox,pos={140,80},size={175,145},title="Misc. Tools"
	Button momentumButton,pos={145,115},size={145,20},title="Convert to momentum",proc=YFSpecButtonProc
	Button expMomButton,pos={145,140},size={145,20},title="Export k Spec",proc=YFSpecButtonProc
	Button expMomButton help={"Copies the k-space map to the current folder"}
	CheckBox checkNaN,pos={160,165},size={75,15},title="0 to NaN",value= 1,proc=YFSpecCheckProc
	CheckBox checkNaN help={"Replaces all zeroes with NaNs in both maps (but not in the cubes themselves)"}

	
	// Make the angle cut display area
	GroupBox angleDispGroupBox,pos={5,270},size={270,270},title="Angle Spec Display"
	Display/W=(10,290,270,535)/HOST=#
	RenameWindow #,GangleMap
	WAVE angleCube = root:Blue:Map:angleCube
	AppendImage angleCube
	ModifyImage angleCube ctab={*,*,ColdWarm,0}
	ModifyGraph wbRGB=(48896,52992,65280), gbRGB=(48896,52992,65280)
	ModifyGraph nticks=4,minor=1,fSize=12,standoff=0,tkLblRot(left)=0,btLen=3
	ModifyGraph margin(left)=30,margin(bottom)=42,margin(top)=13,margin(right)=16
	ModifyGraph swapXY=1
	Label left "\eEnergy (eV)"
	Label bottom "\etheta (deg)"
	SetActiveSubwindow ##
	
	// Make the momentum cut display area
	GroupBox momentumMapDispGroupBox,pos={350,5},size={470,535},title="Momentum Spec Display"
	Display/W=(355,25,805,535)/HOST=# 
	RenameWindow #,GkMap
	WAVE kcube = root:Blue:Map:kcube
	AppendImage kcube
	ModifyImage kcube ctab={*,*,ColdWarm,0}
	ModifyGraph wbRGB=(48896,52992,65280), gbRGB=(48896,52992,65280)
	ModifyGraph swapXY=1
	//ModifyGraph nticks=4,minor=1,fSize=12,standoff=0,tkLblRot(left)=0,btLen=3
	//ModifyGraph margin(left)=50,margin(bottom)=35,margin(top)=5,margin(right)=5
	//ModifyGraph width=0, height={Aspect,DimDelta(kMap,1)*(DimSize(kMap,1)-1)/(DimDelta(kMap,0)*(DimSize(kMap,0)-1))}
	Label left "Energy (eV)"
	Label bottom "k\\Bx\\M (Å\\S-1\\M)"
	SetActiveSubwindow ##
	
	// Get ready to do stuff and update the display
	YFSpecUpdateDisplay()
	
End

Function YFSpecInit()
	
	// Store the current folder location
	String currDF = GetDataFolder(1)
	SetDataFolder root:Blue:Map
	
	// Declare global variables, strings, and waves
	Make /O /T /N=0 listWave
	Make /d/O /N=(2,2) angleCube
	Make /d/O /N=(2,2) kCube
	Variable /G selectedRow = -1
	Variable /G theta
	Variable /G phid
	Variable /G EF=16.745
	Variable /G bin=1
	Variable /G status=0
	Variable /G doNaN = 1
	
	// Return to the current folder
	SetDataFolder currDF
	
	// Get the list of 2D text waves in the folder
	WAVE /T listWave = root:Blue:Map:listWave
	String list = WaveList("*",";","DIMS:2")
	Redimension /N=(ItemsInList(list)) listWave
	listWave = StringFromList(p,list)
	
End

Function YFSpecUpdateDisplay()
	
	NVAR status = root:Blue:Map:status
	
	if (status == 0)	// No DA30 data has been chosen, so disable all the controls
		// Disable the angle cube controls
		SetVariable thetaSetVar disable=2
		SetVariable phidSetVar disable=2
		SetVariable EFSetVar disable=2
		SetVariable binSetVar disable=2
		// Enable some of the misc. tools
		CheckBox checkNaN disable = 2
		Button momentumButton disable = 2
		Button expMomButton disable = 2
		
		SetWindow panelBlueMap#GangleMap, hide=1
		SetWindow panelBlueMap#GkMap, hide=1
		
	elseif (status == 1)	// A momentum map has been created, so display it and enable all the controls
		// Enable the angle cube controls
		SetVariable thetaSetVar disable=0
		SetVariable phidSetVar disable=0
		SetVariable EFSetVar disable=0
		SetVariable binSetVar disable=0
		// Enable the misc. tools
		CheckBox checkNaN disable = 0
		Button momentumButton disable = 0
		Button expMomButton disable = 0
		
		SetWindow panelBlueMap#GangleMap, hide=0
		SetWindow panelBlueMap#GkMap, hide=0
	endif
	
End

Function YFSpecKCube()
   	WAVE angleCube = root:Blue:Map:angleCube
	WAVE kCube = root:Blue:Map:kCube
	NVAR doNaN = root:Blue:Map:doNaN
	NVAR status = root:Blue:Map:status
	NVAR theta = root:Blue:Map:theta
	NVAR phid = root:Blue:Map:phid
	NVAR bin = root:Blue:Map:bin
	NVAR EF = root:Blue:Map:EF
	
   MakeListSpec(angleCube)
   wave momentumlist
   Make /d/O /N=(DimSize(momentumList,0)) tem=0
   tem=momentumlist[p][0]
	// Now find the k-space bounding of the cut
	Variable  kyMin, kyMax
	WaveStats /Q tem
	kyMin = V_min
	kyMax = V_max
	
	variable mapsizey=bin*dimsize(anglecube,1)
	variable delta=(kymax-kymin)/(mapsizey-1)
	
	// Set up the k cube and the Jacobian
	Redimension /N=(DimSize(anglecube,0),mapSizeY) kCube
	SetScale /P y kyMin, delta, kCube
	SetScale /P x DimOffset(anglecube,0), DimDelta(anglecube,0), kCube
	kCube = 0
	duplicate/o kcube, BlueJacobian
	
	// Make the k cube
	Variable pee, i,j
   for(i=0;i<DimSize(anglecube,0);i+=1)
		for (j=i*dimsize(anglecube,1); j<(i+1)*dimsize(anglecube,1); j+=1)
			pee = round((momentumList[j][0] - DimOffset(kCube,1))/DimDelta(kCube,1))
			if (pee >= 0 && pee < DimSize(kCube,1))
				kCube[i][pee] += momentumList[j][1]
				BlueJacobian[i][pee] += 1
			endif
		endfor
		
	endfor
	for(i=0;i<DimSize(kcube,0);i+=1)
	    for(j=0;j<dimsize(kcube,1);j+=1)
	        kCube[i][j] /= max(1,BlueJacobian[i][j])
	    endfor
	endfor
	// Clean up
	KillWaves momentumList, BlueJacobian,tem
	if (doNan)
		kCube =  (kCube[p][q]==0) ? NaN : kCube[p][q]
	endif
	setscale/p x,dimoffset(kcube,0)-EF,dimdelta(kcube,0),kcube
End


function MakeListSpec(angleMap)
  	
	WAVE angleMap
	NVAR theta = root:Blue:Map:theta
	NVAR phid = root:Blue:Map:phid
	
	Make /d/O /N=(DimSize(angleMap,0)*DimSize(angleMap,1),2) momentumList=0

   variable i,j,n=0
   variable ky
   variable phi=0,energy,thetaa
	// Make a list and convert to momentum
	for (i=0; i<DimSize(angleMap,0); i+=1)
	   energy=Dimoffset(angleMap,0)+i*Dimdelta(angleMap,0)
	   
		for (j=0; j<DimSize(angleMap,1); j+=1)
		   
		   thetaa=Dimoffset(angleMap,1)+j*Dimdelta(angleMap,1)
		   
		   ky=0.512*sqrt(energy)*(cos(theta*pi/180)*sin(thetaa*pi/180)-sin(theta*pi/180)*cos((phi+phid)*pi/180)*cos(thetaa*pi/180))

			momentumList[n][0] = ky
          momentumList[n][1] = anglemap[i][j]
			n += 1
		endfor
	endfor
end


Function YFSpecButtonProc(ctrlName)
	String ctrlName
	
	WAVE /T listWave = root:Blue:Map:listWave
	NVAR stage = root:Blue:Map:status
	NVAR doNan = root:Blue:Map:doNan
	NVAR selectedRow = root:Blue:Map:selectedRow
	
	strswitch(ctrlName)
		case "updateButton":
			// Update the wave list
			string list = WaveList("*",";","DIMS:2")
			Redimension /N=(ItemsInList(list)) listWave
			listWave = StringFromList(p,list)
			
			// Reset the status and update the display
			selectedRow = -1
			ListBox specListBox,selRow = selectedRow
			stage = 0
			donan=1
			CheckBox checkNan value = 1
			YFSpecUpdateDisplay()
			break
		case "momentumButton":
			YFSpecKCube()
			YFSpecUpdateDisplay()
			break
		case "expMomButton":
			Duplicate /O root:Blue:Map:kcube $listwave[selectedRow]+"_k"
			break
			Variable /G status = stage
			break
	endswitch
End

Function YFSpecListBoxProc(ctrlName,row,col,event) : ListBoxControl
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	NVAR status = root:Blue:Map:status
	NVAR EF = root:Blue:Map:EF
	WAVE /T specList = root:Blue:Map:specList
	WAVE  anglecube = root:Blue:Map:anglecube
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
	       duplicate/o $listwave[selectedRow],	tempspec
	       redimension/n=(tempspec[0],tempspec[1])  root:Blue:Map:anglecube
	       duplicate/o tempspec, root:Blue:Map:anglecube
	       killwaves tempspec
			YFSpecUpdateDisplay()
		endif
		
	endif
	return 0
End

Function YFSpecCheckProc(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked

	
	strswitch (ctrlName)
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
	endswitch
	
End
