/**
* @author rccomp / http://pericosm.com/
*/

var Tack = (function () {

var thrLc,thrCl //three arrays
var tkLc,tkVl //tack loc,vel Float32Array(*3
var tkCl //tack color        Float32Array(
var tkMa //tack mass         Float32Array(
var tkGp //tack group        Uint16Array(
var tkTy //tack type         Uint16Array(
var tkTm //tack time         Uint16Array(
var tkQa  //tack vector array size
var tkQi  //tack size

var mingrpI, mingrpFc, mingrpLH, mingrpCm
var midgrpI, midgrpFc, midgrpLH, midgrpCm
var medgrpI, medgrpFc, medgrpLH, medgrpCm
var maxgrpI, maxgrpFc, maxgrpLH, maxgrpCm

var mxgrpsz=24, grpszfac=0.75

var gdivs  //grid division vector
var	gdivm  //grid div measure vector
var	gnlow  //grid epsiloned lowbounds vector

var maxgdiv = 256    //max subdivision of space per iteration
var maxtks = 64*1024 //max tacks in simulation
	
function setGroupWorks()
{   	
	var sc = maxtks/4 //16k
	
	//min group - each group of num sc, 8 prt indices
	//16kA    16*9+32*10, 58 bytes per A, total 928 kB
	mingrpI = 
	[ new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ]          
	
	mingrpFc = new Uint16Array( sc ) //status fillcount
	
	mingrpLH = //lowhighbnd 96k words 384k
	[ 
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc )	 ]
	
	mingrpCm = //centerofmass64k words 256k
	[ new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc )  ]
	
	sc = (maxtks/24)|0
	
	//mid group - each group of num sc, 8 prt indices
	midgrpI = //128k words   256k
	[ new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ] 
	
	midgrpFc = new Uint16Array( sc ) //status fillcount
	
	midgrpLH = //lowhighbnd 96k words 384k
	[ 
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc )	 ]
		
	midgrpCm = //centerofmass64k words 256k
	[ new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc )  ]
		
	sc = (maxtks/144)|0
	
	//min group - each group of num sc, 8 prt indices
	medgrpI = //128k words   256k
	[ new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ] 
		
	medgrpFc = new Uint16Array( sc ) //status fillcount
	
	medgrpLH = //lowhighbnd 96k words 384k
	[ 
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ]
		
	medgrpCm = //centerofmass64k words 256k
	[ new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ]
		
	sc = (maxtks/512)|0
	
	//min group - each group of num sc, 8 prt indices
	maxgrpI = //128k words   256k
	[ new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ,
		new Uint16Array( sc ) ] 
		
	maxgrpFc = new Uint16Array( sc ) //status fillcount
	
	maxgrpLH = //lowhighbnd 96k words 384k
	[ 
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ]
		
	maxgrpCm = //centerofmass64k words 256k
	[ new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ,   
		new Float32Array( sc ) ]
		
	//----------------------------------------		
}

//Tack by Sector detail object
var Rsdo = function()
{ return { 
		//vals: new Uint16Array( maxtks ),  //max tacks
		sctr: new Uint16Array( maxgdiv ),   //roster elements
		scts: 0,                        //number of elements
		divs: new Float32Array(3),      //divisions (each axis)
		divm: new Float32Array(3),      //division measure
		lobnd: new Float32Array(3)      //low bound
	}
}

var Rsd=[]; //Tacks by Sector detail array 

function subrostlo3(sec,dvs,dvm,low) //checked
{ var xsec,ysec,zsec

	xsec=sec%dvs[0]
	ysec=((sec/dvs[0])|0)%dvs[1]
	zsec=(sec/(dvs[0]*dvs[1]))|0

  return [low[0]+xsec*dvm[0],low[1]+ysec*dvm[1],low[2]+zsec*dvm[2]]
}

function chopWorld()
{ 
	mkTopRost()
	
	mksubRost( 0, 0, ((mrandom()*(maxgdiv-70)/2)|0)+69 ) 
	
	for( var si=0; si<Rsd[1].scts; si++ )
	{ 
		var tq=Rsd[1].sctr[si+1] - Rsd[1].sctr[si]
		
		if(tq>mxgrpsz)
		{ dosubRost( 1, si, tq) }
		else
		{ doGroup( Rsd[1].sctr[si], Rsd[1].sctr[si+1] )	}
	}
}

//split = tts/tgs
//for uneven spread half groups content > tgs,
//  half groups < tgs
//spread evenness, not efficiently estimable
//so split for smaller group content target
//so morethan half groups < actual maximum group size
//scts = tl/(mxgrpsz*grpszfac)

function dosubRost( flv, fsc, ftq )
{ 
	var mxdivs= (ftq/(mxgrpsz*grpszfac)|0)+1
	mksubRost( flv, fsc, mxdivs ) 
	var lv=flv+1
	
	for(var si=0; si<Rsd[lv].scts; si++)
	{
		var tq=Rsd[lv].sctr[si+1] - Rsd[lv].sctr[si] //tack quantity
		
		if(tq==0){ continue }
		if(tq>mxgrpsz)
		{ dosubRost( lv, si, tq ) }
		else
		{ 
			if( ((si+1)%Rsd[lv].divs[0]!=0) && //next is not in new xline
					( Rsd[lv].sctr[si+2] - Rsd[lv].sctr[si] < mxgrpsz ) 
				) 
			{ doGroup( Rsd[lv].sctr[si], Rsd[lv].sctr[si+2] )
				si++ }
			else
			{ doGroup( Rsd[lv].sctr[si], Rsd[lv].sctr[si+1] ) }
		}
	}
}

function mkTopRost()
{ 
	Rsd[0]=Rsdo()
	Rsd[0].scts=1
	
	Rsd[0].sctr[0]=0                //first tack ix
	Rsd[0].sctr[1]=tkLc.length/3    //after tack ix
	
	var lowx=tkLc[0], lowy=tkLc[1], lowz=tkLc[2]
	var higx=tkLc[0], higy=tkLc[1], higz=tkLc[2]  
	var I=-1;
	
	for(var i=0; i<Rsd[0].sctr[1]; i++)
	{ 
		Rst[i]=i;
		
		if(tkLc[++I]<lowx)    { lowx=tkLc[I] }
		else if(tkLc[I]>higx) { higx=tkLc[I] }
		if(tkLc[++I]<lowy)    { lowy=tkLc[I] } 
		else if(tkLc[I]>higy) { higy=tkLc[I] }
		if(tkLc[++I]<lowz)    { lowz=tkLc[I] } 
		else if(tkLc[I]>higz) { higz=tkLc[I] } 
	}
	
	Rsd[0].lobnd[0]=lowx ; Rsd[0].lobnd[1]=lowy	; Rsd[0].lobnd[2]=lowz	
	Rsd[0].divs[0]=1     ; Rsd[0].divs[1]=1     ; Rsd[0].divs[2]=1
	Rsd[0].divm[0]=higx-lowx 
	Rsd[0].divm[1]=higy-lowy	
	Rsd[0].divm[2]=higz-lowz 


  //hack dither
	var q=1.1
	q=mrandom()*Rsd[0].divm[0]/4
	Rsd[0].divm[0]+=q
	q=mrandom()*Rsd[0].divm[1]/4
	Rsd[0].divm[1]+=q
	q=mrandom()*Rsd[0].divm[2]/4
	Rsd[0].divm[2]+=q
	q=mrandom()*Rsd[0].divm[0]/4
	Rsd[0].lobnd[0]-=q
	q=mrandom()*Rsd[0].divm[1]/4
	Rsd[0].lobnd[1]-=q
	q=mrandom()*Rsd[0].divm[2]/4
	Rsd[0].lobnd[2]-=q
	
	//console.log("Rst was %o", Rst); 

	return
}

var sctsize = new Uint16Array(maxgdiv); //contains quantity of peas in sector
var sctstrt = new Uint16Array(maxgdiv); //contains start index of peas sector in roster
var sctfill = new Uint16Array(maxgdiv); //counts fill of pot while filling roster
var itksec  //contains sector of tack at [ri]
var itkdix  //contains direct index of tack at [ri]
var Rst   //tack indexs by sector

//this functions recursion level involves the density of tacks
//grouping only knows what level to..
//there should not be outliers, the algorithms scale
//begins at level 0

//a roster element contains the start index of a sector in the tacklist
//the tacklist is a list of tack addresses ordered by sector occupation
var dbgc=0,gd=0,bd=0

function mksubRost( flv, fsc, scts ) //lvl uprost index , scts of divs
{ 
	var st = Rsd[flv].sctr[fsc]   //starting in rostertack list
	var ov = Rsd[flv].sctr[fsc+1] //finished in rostertack list
	var lv=flv+1
	
	if(Rsd[lv]==null){ Rsd[lv]=Rsdo() }
		
	var lw3 =subrostlo3(fsc, Rsd[flv].divs, Rsd[flv].divm, Rsd[flv].lobnd)

	calcBestGrid(Rsd[flv].divm, scts)	
	///gdivm, gdivs updated - do not alter 

  //deal with fp precison gunk
	var boge = ( gdivm[0]+gdivm[1]+gdivm[2]  //bogoepsilon
	  +Math.abs(lw3[0])+Math.abs(lw3[1])+Math.abs(lw3[2]) )*Math.sqrt(scts)
	boge=boge/100000000*boge +1/10000000
	//these calibrated to work from tiny to 1000s - check largest scale ...
	 
	
  Rsd[lv].divm[0]=gdivm[0]+boge 
  Rsd[lv].divm[1]=gdivm[1]+boge 
  Rsd[lv].divm[2]=gdivm[2]+boge
	boge=boge/8
	Rsd[lv].lobnd[0]=lw3[0]-boge
	Rsd[lv].lobnd[1]=lw3[1]-boge
	Rsd[lv].lobnd[2]=lw3[2]-boge

  Rsd[lv].divs[0]=gdivs[0]
  Rsd[lv].divs[1]=gdivs[1]
  Rsd[lv].divs[2]=gdivs[2]
	
	Rsd[lv].scts=gdivs[0]*gdivs[1]*gdivs[2]
	scts=Rsd[lv].scts
  
	for(var i=0; i<=scts; i++) { sctfill[i]=0; sctsize[i]=0; }
	var tka, sec, schack=scts*1000
	
	var hib=[  /// /// ///
	  lw3[0]+Rsd[flv].divm[0],
	  lw3[1]+Rsd[flv].divm[1],
	  lw3[2]+Rsd[flv].divm[2]
		]  /// /// ///
		
	for(var ri=st; ri<ov; ri++)  //for through stretch to make note of members sectors
	{ tka=Rst[ri]*3
		sec=locToSector(
		     tkLc[tka], tkLc[tka+1], tkLc[tka+2],
				 Rsd[lv].lobnd, Rsd[lv].divm, Rsd[lv].divs );
	  
    sec=(sec+schack)%scts		
		sctsize[sec]++;
		
		if(!(sec>-1&&sec<scts))
		//~ if(0)
		{ bd++;
		  console.log( "g=%s b=%s sc=%s fl=%s xl=%s xi=%s xh=%s yl=%s yi=%s yh=%s zl=%s zi=%s zh=%s",gd,bd,sec,flv, 
			lw3[0],tkLc[tka],hib[0],lw3[1],tkLc[tka+1], hib[1],lw3[2], tkLc[tka+2], hib[2] 
			);
		}
	  else{ gd++ }
		
		itksec[ri]=sec;      //the sector of the tack at the tacksector
		itkdix[ri]=Rst[ri];  //the direct index of the tack at roster indirect
		
	}
	
	//each of st to ov in Rst, sector is noted in itksec, tacki is note in itkdix
	
	var sst=Rsd[flv].sctr[fsc];  //sector start in Rst
	for(var sec=0; sec<=scts; sec++)  //loop through sector to make roster index
	{ Rsd[lv].sctr[sec]=sst; sst+=sctsize[sec]; }
	
	for(var ri=st;ri<ov;ri++)  //loop through tack refs for put into roster
	{ sec=itksec[ri]  
		Rst[ Rsd[lv].sctr[sec]+(sctfill[sec]++) ]
		  =itkdix[ri]; 
	}
	
	return
}

function locToSector(x,y,z,lw,dvm,dvs)
{ return (
    (((x-lw[0])/dvm[0])|0)                      //*1
	+ ((((y-lw[1])/dvm[1])|0)*dvs[0])             //*1*x 
	+ ((((z-lw[2])/dvm[2])|0)*dvs[0]*dvs[1])      //*1*x*y 
	) 
}  

//   ((((((z-lw[2])/dvm[2])|0)*dvs[1])  //alt order
//+  (((y-lw[1])/dvm[1])|0))*dvs[0])
//+  (((x-lw[0])/dvm[0])|0)

var fsz3=0,fmxdiv=0

function calcBestGrid(sz3,mxdiv)
{ 	
  if((mxdiv==fmxdiv)&&fsz3[0]==sz3[0]  //simple caching last calc
	  &&fsz3[0]==sz3[0]&&fsz3[0]==sz3[0]) { return }
	var put3=[sz3[0],sz3[1],sz3[2]]
	//put3[0]=put3[0]*1.6  //-squeeze more into this axis
	put3=sadd3( put3, (put3[0]+put3[1]+put3[2]+0.1)/1000.0)
	var vol=put3[0]*put3[1]*put3[2]; 
	var fak= Math.pow( mxdiv /vol , (0.333334) ); //multi by cbrt to make vol mxdiv
	put3=smult3(put3,fak)            //scalar mult by scale to mxdiv volume
	var rnk=sorti012(put3);
	
	//boost the lowest value up if low, avoid error
	if((put3[rnk[0]])<(mxdiv/65|0)+1) {	 
		var lv=((put3[rnk[0]]|0)+1)
		fak=Math.sqrt(put3[rnk[0]]/lv)
		put3[rnk[0]]=lv
		put3[rnk[1]]=put3[rnk[1]]*fak
		put3[rnk[2]]=put3[rnk[2]]*fak
		if (put3[rnk[1]]<put3[rnk[0]]) { //second has become lower 
			fak=put3[rnk[1]]/(put3[rnk[0]])
			put3[rnk[1]]=put3[rnk[0]] 
			put3[rnk[2]]=put3[rnk[2]]*fak 
		} 
	}
	put3=[ put3[0]|0, put3[1]|0, put3[2]|0 ]
	fak=Math.sqrt((mxdiv/put3[0]*put3[1]*put3[2]))
	if(fak<0.95)
	{ if( put3[rnk[0]]*((put3[rnk[1]]*fak-1)|0)*((put3[rnk[2]]*fak)|0)<=mxdiv )
		{  put3[rnk[1]]=((put3[rnk[1]]*fak-1)|0)
			put3[rnk[2]]=((put3[rnk[2]]*fak)|0) } }
	if( put3[rnk[0]]*((put3[rnk[1]]+1)|0)*((put3[rnk[2]])|0)<=mxdiv )
	{ put3[rnk[1]]+=1; } 
	
	put3[rnk[2]]=(mxdiv/(put3[rnk[0]]*put3[rnk[1]]))|0;
	
	if(put3[rnk[2]]<put3[rnk[1]])
	{ fak=put3[rnk[2]]; put3[rnk[2]]=put3[rnk[1]]; put3[rnk[1]]=fak }
	
	//print(put3, dlen3 , put3[0], put3[0]*put3[1])
	//put3 is now split numbers of grid
	gdivm =[ sz3[0]/put3[0], sz3[1]/put3[1], sz3[2]/put3[2] ];
	gdivs =put3

	return 
}

function sorti012(m) 
{ if(m[0]<m[1])                               
	{ if(m[1]<m[2]){ return [0,1,2]; }  //1 not largest  
		if(m[0]<m[2]){ return [0,2,1]; }
		return [2,0,1]  }                                     
	if(m[0]<m[2]) { return [1,0,2]; } 
	if(m[1]<m[2]) { return [1,2,0]; }
	return [2,1,0]; 
}

function delt3(a,b)
{ return [b[0]-a[0], b[1]-a[1], b[2]-a[2]] }
function smult3(a,b)
{ return [ a[0]*b, a[1]*b, a[2]*b] }
function sadd3(a,b)
{ return [ a[0]+b, a[1]+b, a[2]+b] }

function delt3r(a,b,r)
{ r[0]=b[0]-a[0];r[1]=b[1]-a[1];r[2]=b[2]-a[2];return r }
function smult3r(a,b,r)
{ r[0]=a[0]*b;r[1]=a[1]*b;r[2]=a[2]*b;return r }
function sadd3r(a,b,r)
{ r[0]=a[0]+b;r[1]=a[1]+b;r[2]=a[2]+b;return r }


//functions which sweep prime addressing space
//functions which react on time integrating - real time


//travelling particles cut across the space which its neibours
//differ by, faster than the spacing
//particles interaction scope is 
//travelling particles extend their groups positional span
//positional span increases likelyness of contact

//splitting cuboid
//given a cuboid volume, knowing density
//calculating string travel in each dimension
//on travel being greater than density*splitfactor
//change group > compare travel to other groups
//and select on better than (accept factor)
//or make new group

//splitfactor depends on desired grpsize, cuboid population, and 
//  decisionscore = (decisionscore +1) /2  marginal contribution
//  decisionscore = (decisionscore +1/2)   material contribution  

//running center  
//  rn_cent+=pdif*(pmass/rn_mass)

//keeping a sum score for each pot
//how far is it from testees? hehe
//how does that compare to global average distances
//how would it differ to accept modified on individual pot scores?
//pot with high score should accept on lesser fit because it is far?
//or greater because it is far?,
//potscoring may be not helpful
//this process is designed to do adequate grouping quickly
//not optimum grouping...

//objects for doGroup
var opos = { x:0, y:0, z:0 }
var hpos = { x:0, y:0, z:0 }
var lpos = { x:0, y:0, z:0 }
var ovel = { x:0, y:0, z:0 }
var opsv = { x:0, y:0, z:0 }       //pos+vel
var potGc= { x:0, y:0, z:0, w:0 }

var potsmax =16,potsize=8
var ptfill = new Array(potsmax)   //pot numbers preGn[0..7]
var pttaks = [ new Uint16Array(potsize), //pot index preGi[0..7]
new Uint16Array(potsize),
new Uint16Array(potsize),
new Uint16Array(potsize)] //...potmax
var ptcntr = [potsmax][6] //pot center pos&vel preGc[g][0..2]	

function doGroupx(rst,rov) //stub to test chopworld
{ var tk,tka,r,g,b
  r=((2+rst)%7)/6 
	g=((1+rst)%6)/5 
	b=((2+rst)%5)/4
	
	for(var ri=rst; ri<rov; ri++)  //for through stretch to make note of members sectors
	{ tk=Rst[ri]; tka=tk*3;
	  tkCl[tka]=r; tkCl[tka+1]=g; tkCl[tka+2]=b;
	  //tkLc[tka]=r; tkLc[tka+1]=g; tkLc[tka+2]=b;
  }	
	return
}



function doGroupy(rst,rov) //stub to test chopworld
{ avgLc[0]=0.0;avgLc[1]=0.0;avgLc[2]=0.0;
	gsz=rov-rst
	
	for(var ri=rst; ri<rov; ri++)  //for through stretch to make note of members sectors
	{ tk=Rst[ri]; tka=tk*3;
	  avgLc[0]+=tkLc[tka++];avgLc[1]+=tkLc[tka++];avgLc[2]+=tkLc[tka]  
  }
	avgLc[0]/=gsz; avgLc[1]/=gsz; avgLc[2]/=gsz
	var ff=100000
	for(var ri=rst; ri<rov; ri++)  //for through stretch to make note of members sectors
	{ tk=Rst[ri]; tka=tk*3;
	  dfLc=tkLc[tka]-avgLc[0]
		if(dfLc>0){ tkVl[tka]+=Math.sqrt(dfLc)/ff }
		else    { tkVl[tka]-=Math.sqrt(-dfLc)/ff }
	  dfLc=tkLc[++tka]-avgLc[1]
		if(dfLc>0){ tkVl[tka]+=Math.sqrt(dfLc)/ff }
		else    { tkVl[tka]-=Math.sqrt(-dfLc)/ff }
	  dfLc=tkLc[++tka]-avgLc[2]
		if(dfLc>0){ tkVl[tka]+=Math.sqrt(dfLc)/ff }
		else    { tkVl[tka]-=Math.sqrt(-dfLc)/ff }
  }	
}

var avgLc =[0.0, 0.0, 0.0]
var avgVl =[0.0, 0.0, 0.0]
var dfLc, gsz, tka,tkb

function doGroupf(rst,rov) //stub to test chopworld
{
  return
}


function doGroup(rst,rov) //stub to test chopworld
{ 
  //~ avgLc[0]=0.0; avgLc[1]=0.0; avgLc[2]=0.0;
  avgVl[0]=0.0; avgVl[1]=0.0; avgVl[2]=0.0;
	gsz=rov-rst
	
	for(var ri=rst; ri<rov; ri++)  //ds avg loc and avg vl
	{ tka=Rst[ri]*3;
	  avgVl[0]+=tkVl[tka]; tka++
		avgVl[1]+=tkVl[tka]; tka++
		avgVl[2]+=tkVl[tka]  
		//~ avgLc[0]+=tkLc[tka];
		//~ avgLc[1]+=tkLc[tka];
		//~ avgLc[2]+=tkLc[tka];
  }	
	
	avgVl[0]/=gsz; avgVl[1]/=gsz; avgVl[2]/=gsz
	//~ avgLc[0]/=gsz; avgLc[1]/=gsz; avgLc[2]/=gsz
	
	//~ for(var ri=rst; ri<rov; ri++)  //
	//~ { tka=Rst[ri]*3;
	  //~ dfLc=(tkLc[tka]-avgLc[0])*(tkLc[tka]-avgLc[0]); tka++
	  //~ dfLc+=(tkLc[tka]-avgLc[1])*(tkLc[tka]-avgLc[1]); tka++
	  //~ dfLc+=(tkLc[tka]-avgLc[2])*(tkLc[tka]-avgLc[2]);
		//~ dfLc=Math.sqrt(dfLc)+0.01; 
		//~ tka-=2;
		//~ var dfLc3=dflc
		//~ // tkVl[tka]=(tkVl[tka]*799 + avgVl[0] )/800 ;tka++
		//~ // tkVl[tka]=(tkVl[tka]*799 + avgVl[1] )/800 ;tka++
		//~ // tkVl[tka]=(tkVl[tka]*799 + avgVl[2] )/800
		//~ tkVl[tka]=(tkVl[tka]*dfLc + avgVl[0] )/(dfLc+1) ;tka++
		//~ tkVl[tka]=(tkVl[tka]*dfLc + avgVl[1] )/(dfLc+1) ;tka++
		//~ tkVl[tka]=(tkVl[tka]*dfLc + avgVl[2] )/(dfLc+1)
  //~ }	
	var dd=0.02
	
	function signt(n){ return n<0 ? -1:1 }
	
	var p
	var rr=rov-rst
	
	for(var ri=rst; ri<rov; ri++)  //
	{ tka=Rst[ri]*3;
	  tkb=Rst[rst+((ri-rst+1)%rr)]*3;
				
	  p=(tkLc[tka]-tkLc[tkb++]) 
		p+=signt(p)  
	  tkLc[tka++]+=dd/p;
		
	  p=(tkLc[tka]-tkLc[tkb++]) 
		p+=signt(p)  
	  tkLc[tka++]+=dd/p;
		
	  p=(tkLc[tka]-tkLc[tkb]) 
		p+=signt(p)  
	  tkLc[tka++]+=dd/p;
		
		tka-=3;
		//~ var dfLc3=dflc
	  tkVl[tka]=(tkVl[tka]*799 + avgVl[0] )/800 ;tka++
	  tkVl[tka]=(tkVl[tka]*799 + avgVl[1] )/800 ;tka++
	  tkVl[tka]=(tkVl[tka]*799 + avgVl[2] )/800
		//tkVl[tka]=(tkVl[tka]*dfLc + avgVl[0] )/(dfLc+1) ;tka++
		//tkVl[tka]=(tkVl[tka]*dfLc + avgVl[1] )/(dfLc+1) ;tka++
		//tkVl[tka]=(tkVl[tka]*dfLc + avgVl[2] )/(dfLc+1)
  }
	
}

function groupRostx(rsbeg,rsfin,roster)
{   
	function pottogrp(grpI,pot)
	{
		var gpin=mingrpFc[grpI]
		var r=0, ix=0, iy=0, iz=0
		
		if(mingrpFc[grpI]==0) //init.reset group
		{ mingrpCm[grpI][0]=0; mingrpCm[grpI][1]=0
			mingrpCm[grpI][2]=0; mingrpCm[grpI][3]=0
		}
		
		for(var p=0; p<ptfill[pot]; p++)
		{ var t=pttaks[p]
			var ix=t*3, iy=t*3+1, iz=t*3+1		
			
			opos = { x:tkLc[ix], y:tkLc[iy], z:tkLc[iz] } //another syntax may be faster
			hpos = { x:tkLc[ix], y:tkLc[iy], z:tkLc[iz] }
			lpos = { x:tkLc[ix], y:tkLc[iy], z:tkLc[iz] }
			ovel = { x:tkVl[ix], y:tkVl[iy], z:tkVl[iz] }						
			opsv = { x:opos.x+ovel.x, y:opos.y+ovel.y, z:opos.z+ovel.z }
			
			r = tkD[t] //radii to calc bounding box
			if(opos.x<opsv.x) { hpos.x=opsv.x+r ; lpos.x=lpos.x-r }
			else{ lpos.x=opos.x+r ; hpos.x=hpos.x-r }
			if(opos.y<opsv.y) { hpos.y=opsv.y+r ; lpos.y=lpos.y-r }
			else{ lpos.y=opos.y+r ; hpos.y=hpos.y-r }
			if(opos.z<opsv.z) { hpos.z=opsv.z+r ; lpos.z=lpos.z-r }
			else{ lpos.z=opos.z+r ; hpos.z=hpos.z-r }
			
			if(lpos.x>mingrpLH[grpI][0]) { lpos.x=mingrpLH[grpI][0] }
			if(lpos.y>mingrpLH[grpI][1]) { lpos.y=mingrpLH[grpI][1] }
			if(lpos.z>mingrpLH[grpI][2]) { lpos.z=mingrpLH[grpI][2] }
			
			if(hpos.x<mingrpLH[grpI][3]) { hpos.x=mingrpLH[grpI][3] }
			if(hpos.y<mingrpLH[grpI][4]) { hpos.y=mingrpLH[grpI][4] }
			if(hpos.z<mingrpLH[grpI][5]) { hpos.z=mingrpLH[grpI][5] }
			
			potGc.x += tkLc[ix] * tkMa[t]
			potGc.y += tkLc[iy] * tkMa[t]
			potGc.z += tkLc[iz] * tkMa[t]
			potGc.w += tkMa[t]          
			
			tkGp[t]=grpI
			mingrpI[grpI][gpin++] = t
		}//finnish taking tacks 
		
		var gw = mingrpCm[grpI][3]  //group weight
		var tw = potGc.w+gw         //total weight
		
		//save grp data 
		mingrpFc[grpI] = gpin
		
		mingrpCm[grpI][0] = (potGc.x*potGc.w + mingrpCm[grpI][0]*gw) /tw 
		mingrpCm[grpI][1] = (potGc.y*potGc.w + mingrpCm[grpI][1]*gw) /tw 
		mingrpCm[grpI][2] = (potGc.z*potGc.w + mingrpCm[grpI][2]*gw) /tw 
		mingrpCm[grpI][3] = tw
		
		mingrpLH[grpI][0] = lpos.x ; mingrpLH[grpI][3] = hpos.x 
		mingrpLH[grpI][1] = lpos.x ; mingrpLH[grpI][4] = hpos.x
		mingrpLH[grpI][2] = lpos.x ; mingrpLH[grpI][5] = hpos.x
		
	}
	
	function poptopot(Pt,Tk)
	{
		if(ptfill[Pt]==4)
		{ Pt=++freepot; 
			if(freepot>16){ laspot=16; } ///clamp dangerous
		}
		pttaks[Pt]=Tk; //pot index preGi[0..7]
		ptfill[Pt]++;  //pot numbers preGn[0..7]
		var q=ptfill[Pt],p=ptfill[Pt]-1 
		Tk=Tk*3
		ptcntr[Pt][0]=(ptcntr[Pt][0]*p+tkLc[Tk]  )/q
		ptcntr[Pt][1]=(ptcntr[Pt][1]*p+tkLc[Tk+1])/q
		ptcntr[Pt][2]=(ptcntr[Pt][2]*p+tkLc[Tk+2])/q
		ptcntr[Pt][3]=(ptcntr[Pt][3]*p+tkVl[Tk]  )/q
		ptcntr[Pt][4]=(ptcntr[Pt][4]*p+tkVl[Tk+1])/q
		ptcntr[Pt][5]=(ptcntr[Pt][5]*p+tkVl[Tk+2])/q
	}
	
	function cal3dist(Tk,px,py,pz,vx,vy,vz)
	{ 
		dist3=
		Math.abs(tkLc[Tk]-px) + Math.abs(tkLc[Tk+1]-py) + Math.abs(tkLc[Tk+2]-pz)
		+(Math.abs(tkLc[Tk  ]+tkVl[Tk]/2 - px+vx/2)  
			+ Math.abs(tkLc[Tk+1]+tkVl[Tk+1]/2 - py+vy/2)
			+ Math.abs(tkLc[Tk+2]+tkVl[Tk+2]/2 - pz+vz/2))*2;
		
		if(d3prev/dist3>1.1)
		{ d3ups+=dist3 ; d3upnum++ //gotten bigger 
		}else{
			d3dns+=dist3 ; d3dnnum++ //gotten smaller
		}
		d3prev=dist3
		return
		
		//purpose of measurement is to favour certain relationships
		//between tacks. nearness/farness, approachingness/departingnesss
		//high speed usually means disconnect unless there is aggreeance
	}
	
	function cal3distpot(pta,ptb)
	{ 
		dist3=
		Math.abs(ptcntr[pta][0]-ptcntr[ptb][0]) 
		+ Math.abs(ptcntr[pta][1]-ptcntr[ptb][1])
		+ Math.abs(ptcntr[pta][2]-ptcntr[ptb][2])
		+(Math.abs(ptcntr[pta][3]-ptcntr[ptb][3])
			+ Math.abs(ptcntr[pta][4]-ptcntr[ptb][4])
			+ Math.abs(ptcntr[pta][5]-ptcntr[ptb][5]))*2
		
		if(d3prev/dist3>1.1)       //no longer useful for pot matching?
		{ d3ups+=dist3 ; d3upnum++ //gotten bigger 
		}else{
			d3dns+=dist3 ; d3dnnum++ //gotten smaller
		}
		d3prev=dist3
		return
		
		//purpose of measurement is to favour certain relationships
		//between tacks. nearness/farness, approachingness/departingnesss
		//high speed usually means disconnect unless there is aggreeance
	}
	
	var rsbeg= rostix[rnum]
	var rsnxt= rsbeg
	var rsfin= rostix[rnum+1]
	var rslen= rsbeg -rsfin //roster index convention
	//densityv[0]=population/divpvect[0]
	
	var d3prev=0
	
	//initialise first tacks to pots, observe distance calcs
	//but kick off d3ups,d3dwns,dist3,d3pre with 
	//estimates from basic density of cuboid
	//ensure roster order is shuffled with suffle sequence
	//the 3dist cals effectively shrink the bounds of the fit
	//the low 3dist cal will still be much bigger than group accept
	//group accept is ~ 3distcal / pop
	//the greater the pop, the smaller the grouping value
	
	var sepb=(bndx+bndy+bndz)
	var lod3=0xffffffff,dntarg,d3dnsum=sepb,d3upsum=sepb,d3upnum=2,d3dnnum=2
	var pots=pttaks.length,freepot=0
	var potTarg=(1+rslen/5)|0
	
	while(rsnxt<potTarg)
	{ 
		var i=roster[ (rsnxt+rslen-1)%rslen ]*3
		
		cal3dist( roster[rsnxt]*3,
			tkLc[i], tkLc[i+1], tkLc[i+2], 
			tkVl[i], tkVl[i+1], tkVl[i+2] )
		
		poptopot(freepot++,roster[rsnxt++]) 
	}
	
	//fill likepots to max depth of 4
	while (rsnxt<rsfin) 
	{ lod3=0xffffffff
		for(var wp =0; wp<freepot; wp++)
		{ if(ptfill[wp]<plen)
			{
				cal3dist(roster[rsnxt]*3,
					ptcntr[wp][0],ptcntr[wp][1],ptcntr[wp][2], 
					ptcntr[wp][3],ptcntr[wp][4],ptcntr[wp][5] )
				
				if( dist3 < lod3 )
				{ lod3=dist3; gopo=wp; }
				
				//d3sum is not required when tracking up and down
				//d3sum, d3num, d3ups, d3dns d3upnum, d3dnnumd
				// (d3ups/d3upnum) / (d3dns/d3dnnum) -- compare to d3sum 
				//the separation measurement of larger pairs and smaller pairs
				//informs of uniformity of spread of points,
				//it settles on a value which is a natural factor 
				//separated by actual diversity of locations
				//that factor should be calibrated by
				//testing on types of uniform spread tacks
				//if self setting, the meaning becomes more,
				//complex and obscure open to new flaws
				//the diff between up and down, informs
				//chancieness of measurement
				//when smaller, measurement is safer,
				//so increase acceptance...
			}
		}
		
		//upsm = upsm + dwsm * (un-dn)/(un+dn)
		var dntarg = (d3dnsum + d3upsum *(d3upnum-d3dnnum)/(d3upnum+d3dnnum))/d3dnnum
		//var uptarg = (d3upsum + d3dnsum *(d3upnum-d3dnnum)/(d3upnum+d3dnnum))/d3upnum
		
		if( lod3 > dntarg*fac2/rslen )
		{ if((++freepot)==pots) freepot-- 
				poptopot(freepot,roster[rsnxt++]); }
		else
		{ d3dnsum-=lod3; d3dnnum-- 
			poptopot(gopo,roster[rsnxt++]); 
		}
		
	}//all tacks in pots
	
	var pota=0
	while (pota<freepot)
	{
		pottogrp(wkgrp,pota)
		
		for(var potb=pota+1;potb<pots;potb++)
		{ if(ptfill[potb]+mingrpFc[wkgrp]<grplen)
			{ cal3distpot(pota,potb)
				
				if(dist3< (fc*(d3upsum+d3dnsum *(d3upnum-d3dnnum)/(d3upnum+d3dnnum))/d3upnum)/rslen )
				{ pottogrp(wkgrp,potb)
					if(mingrpFc[wkgrp]>=grplen){ wkgrp++; break }			
				}
			}
		}
	}
}


function initThrees(loca,colr)
{ thrLc=loca; thrCl=colr; }

function syncTacks()
{
  for(var i=0;i<thrLc.length;i++)
	{ thrLc[i]=tkLc[i] }
  for(var i=0;i<thrCl.length;i++)
	{ thrCl[i]=tkCl[i]; }
}

function initTacks()
{
	tkQa=thrLc.length
	tkQi=tkQa/3
	
	tkLc = new Float32Array(tkQa)
	tkVl = new Float32Array(tkQa)
	tkCl = new Float32Array(tkQa)
	tkMa = new Float32Array(tkQi)
	tkGp = new Uint16Array(tkQi)
	tkTy = new Uint16Array(tkQi)
	tkTm = new Uint16Array(tkQi)
	
	itksec  = new Uint8Array (tkQi); //contains sector of tack at [ri]
	itkdix  = new Uint16Array (tkQi); //contains address of tack at [ri]
	Rst = new Uint16Array(tkQi); //contains tack

	initPhys2(tkLc,90)
  initPhys(tkVl,15)

  var r=0,g=0,b=0;
  for(var i=0;i<tkQa; )
	{ 
	  tkCl[i++]=r/1000; r=(r+20)%1000 
	  tkCl[i++]=g/1000; g=(g+100)%1000 
	  tkCl[i++]=b/1000; b=(b+950)%1000 
	}

	function initPhys2(phys,spread)
	{
		var i=0
		
		while(i < phys.length)
		{ 

			phys[i++]=(Math.random()-0.5)*spread
			phys[i++]=(Math.random()-0.5)*spread
			phys[i++]=(Math.random()-0.5)*spread
			
		}
	}

	function initPhys(phys,spread)
	{
		var i=0, n=spread, n2=n/2; 
		var bnc=0, bncp=0, bncpx=0, bncpy=0, bncpz=0
		var x,y,z, wx=0,wy=0,wz=0
		
		while(i < phys.length)
		{ if(bnc<3)
			{ 
				bnc=(mrandom()*mrandom()*(phys.length-i)/5)+1
				bncpx=(((mrandom()-0.5)*n)*3)/4
				bncpy=(((mrandom()-0.5)*n)*3)/4
				bncpz=(((mrandom()-0.5)*n)*3)/4
				wx=mrandom()*n2
				wy=mrandom()*n2
				wz=mrandom()*n2
			}
			
			x=(mrandom()-0.5)*wx; y=(mrandom()-0.5)*wy; z=(mrandom()-0.5)*wz

			phys[i]=bncpx+x; phys[i+1]=bncpy+y; phys[i+2]=bncpz+z

			i=i+3
			bnc--
		}
	}
}

var shrinki=0;
function shrink(t)
{
  shrinki++;
	var shfc= (1+Math.sin(shrinki/100))*t
	for(var i=0;i<tkQa;)
	{ 
		tkVl[i]-=tkLc[i++]*shfc
	}
}


var shrinki2=0;
function shrink2(f)
{
  shrinki2++;
	var shfc= (Math.sin(shrinki2/88))*f
	for(var i=0;i<tkQa;)
	{ 
		tkLc[i]=tkLc[i++]*shfc
		tkLc[i]=tkLc[i++]*shfc
		tkLc[i]=tkLc[i++]*shfc
	}
}


function velmove(f)
{
  for(var i=0;i<tkQa;)
	{ 
		tkLc[i]=tkLc[i]+tkVl[i++]*f
		tkLc[i]=tkLc[i]+tkVl[i++]*f
		tkLc[i]=tkLc[i]+tkVl[i++]*f
	}
}

function velcolor(f)
{
  for(var i=0;i<tkQa;)
	{ 
		tkCl[i]=0.15+0.85/(Math.pow(Math.max(1,Math.abs(tkVl[i++])),f) )
		tkCl[i]=0.15+0.85/(Math.pow(Math.max(1,Math.abs(tkVl[i++])),f) )
		tkCl[i]=0.15+0.85/(Math.pow(Math.max(1,Math.abs(tkVl[i++])),f) )
	}
}
	
	//move all tacks by their velocity
  //have an elastic central tether, 
	//whose force increases by square of tack distance from center
	//this will arrest all tacks at points where force balances out
	//modulate the power of the elestic force to keep tacks in motion



function rndmovep()
{
	var scale=1;
	for ( var i = 0; i < rc.length; i += 3 ) {
		// rawpos
		var vx = mrandom();
		var vy = mrandom();
		var vz = mrandom();
		tkVl[ i ]     = vx;
		tkVl[ i + 1 ] = vy;
		tkVl[ i + 2 ] = vz;
		tkLc[ i ]     += (vx*scale-scale/2)|0;
		tkLc[ i + 1 ] += (vy*scale-scale/2)|0;
		tkLc[ i + 2 ] += (vz*scale-scale/2)|0;
	}
	return
}

function blasta()
{ 
  //console.log("$s ",mrandom())
  return }

var qRnd=0xcb8c29e3
function mrandom() ///shr3 32
{ qRnd ^= (qRnd << 17); qRnd ^= (qRnd >>> 13); qRnd ^= (qRnd << 5); 
	return (qRnd&0x7fffffff)/0x7fffffff }
		
return {
	initThrees     :  initThrees,
	initTacks      :  initTacks,
	syncTacks      :  syncTacks,
	shrink         :  shrink,
	shrink2        :  shrink2,
	velmove        :  velmove,
	velcolor       :  velcolor,
	chopWorld      :  chopWorld,
	blasta         :  blasta,
	setGroupWorks  :  setGroupWorks
};})();


//~ function sumdiff(vert0,vert1)
//~ {}
//~ function rdvert(varry, offset)  //read vert
//~ { return varry[offset]; }
//~ function rwvert(varry, offset, val)  //write vert
//~ { varry[offset]=val; }
//~ function deltvert(varry, offset0, offset1)
//~ { return [a,b,c] }
//~ function delt(varry, off0,off1)
//~ { return varry[off0]-varry[off0]; }


