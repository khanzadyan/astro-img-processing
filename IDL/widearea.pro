pro widearea

	; refimg=readfits('./FIR/iras100map.fits',refimgH)
	refimg=readfits('./NIR/NB2122/nb2122wcs.fits',refimgH)

	iras4img=readfits('./FIR/iras100map.fits',iras4imgH)
	
	avmapimginit=readfits('./NIR/avmap49_-10to0.fits',avmapimginitH)
	; avmapimginit=readfits('./NIR/av_jhhk_4arcsec_paper.fits',avmapimginitH)
	; avmapimginit=readfits('./Optical/78023543.fits',avmapimginitH)
	; optimg=readfits('./Optical/dss2-red.fits',optimgH)
	nirimg=readfits('./NIR/NB2122/nb2122wcs.fits',nirimgH)
	
	comapimg=readfits('./Radio/co-region.fits',comapimgH)
	
	avmapimg=avmapimginit(356:650,326:619)
	hextract,avmapimginit,avmapimginitH,avmapimg,avmapimgH,356,650,326,619
	
	; avmapimg=avmapimginit(1723:1803,798:898)
	; hextract,avmapimginit,avmapimginitH,avmapimg,avmapimgH,1723,1803,798,898
	heuler,avmapimgH,/CELESTIAL
	
	
	; hastrom,optimg,optimgH,refimgH,MISSING=0
	hastrom,nirimg,nirimgH,refimgH,MISSING=0
	hastrom,avmapimg,avmapimgH,refimgH,MISSING=0
	hastrom,iras4img,iras4imgH,refimgH,MISSING=0
	hastrom,comapimg,comapimgH,refimgH,MISSING=0
	
	; 05 40 24.233 +23 50 54.62
	; obj_ra=15*ten(05,40,24.233)
	; obj_dec=ten(23,50,54.62)

	obj_ra=15*ten(05,40,29.994)
	obj_dec=ten(23,52,12.81)


	
	; defining the center of the field
	adxy,refimgH,obj_ra,obj_dec,obj_x,obj_y
	HalfLength=900
	Start_x=obj_x - HalfLength
	End_x=obj_x + HalfLength
	Start_y=obj_y - HalfLength
	End_y=obj_y + HalfLength
	
	
	; opt=optimg(Start_x:End_x,Start_y:End_y)
	; hextract,optimg,optimgH,opt,optH,Start_x,End_x,Start_y,End_y
	; sky,opt,opt_skyvalue,opt_skysig,/NAN
	
	nir=nirimg(Start_x:End_x,Start_y:End_y)
	hextract,nirimg,nirimgH,nir,nirH,Start_x,End_x,Start_y,End_y
	sky,nir,nir_skyvalue,nir_skysig,/NAN
	
	iras=iras4img(Start_x:End_x,Start_y:End_y)
	hextract,iras4img,iras4imgH,iras,irasH,Start_x,End_x,Start_y,End_y
	sky,iras,iras_skyvalue,iras_skysig,/NAN
	
	avmapimg=smooth(avmapimg,200,/NAN)
	avmap=avmapimg(Start_x:End_x,Start_y:End_y)
	hextract,avmapimg,avmapimgH,avmap,avmapH,Start_x,End_x,Start_y,End_y
	sky,avmap,avmap_skyvalue,avmap_skysig,/NAN
	; avmap=smooth(avmap,200,/NAN)
	
	comapimg=smooth(comapimg,200,/NAN)
	comap=comapimg(Start_x:End_x,Start_y:End_y)
	hextract,comapimg,comapimgH,comap,comapH,Start_x,End_x,Start_y,End_y
	sky,comap,comap_skyvalue,comap_skysig,/NAN
	sky,comap,comap_skyvalue,comap_skysig,/NAN
	; comap=smooth(comap,3,/NAN)

	naxis1 = sxpar(irasH, 'NAXIS1')
	naxis2 = sxpar(irasH, 'NAXIS2')


; Reading the catalougs ==================================================

	FMTSpitzer='A,I,I,F,I,I,F,I,F,F'
	extast,nirH,nir_astrometry ; extracting the astrometry info
	
	READCOL,'catalogs/spitzer-complete.txt',F=FMTSpitzer,$
		sp_name,sp_hr,sp_min,sp_sec,sp_deg,sp_dmin,sp_dsec,sp_class,sp_kv,sp_alpha
	print,'reading spitzer objects....'

	spcount=size(sp_name)
	last_sp=spcount(1)
	sp_ra=dindgen(last_sp)
	sp_dec=dindgen(last_sp)
	sp_x=fltarr(last_sp)
	sp_y=fltarr(last_sp)

	for p=0, last_sp-1 do begin
		sp_ra[p]=15*ten(sp_hr[p],sp_min[p],sp_sec[p])
		sp_dec[p]=ten(sp_deg[p],sp_dmin[p],sp_dsec[p])
	endfor

	ad2xy, sp_ra, sp_dec, nir_astrometry, sp_x, sp_y
	; writefits,'avmap.fits',avmap,avmapH

; ======================== Ploting parameters ===================

	!X.MARGIN=0
	!Y.MARGIN=0
	!X.THICK=3
	!Y.THICK=3
	!P.CHARSIZE=-1
	; 
	; sizeiras=size(iras)
	; sizex=sizeiras(1)
	; sizey=sizeiras(2)
		
	cgLoadct,0, /Reverse
	; !P.MULTI=[0,1,1]
	!P.POSITION=[0.11,0.09,0.96,0.94]
	!p.font=0 ; nice post script fonts
	set_plot,'ps',/copy, /interpolate
	
	device,filename='widearea.eps',/encapsulated,/color,bits=8, $
	/helvetica,/isolatin1,$
	xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	; axes=lonarr(3)
	
	; final_image=lonarr(3,sizex,sizey)
	; ; construction of rgb image from the suplied min and max
	; final_Image[0,*,*] = bytscl(comap,min=(comap_skysig),max=max(comap))
	; final_Image[1,*,*] = bytscl(iras,min=(iras_skysig*5),max=(iras_skysig*12))
	; final_Image[2,*,*] = bytscl(avmap,min=(avmap_skysig),max=max(avmap))
	; 
	; tvscale,final_image,position=!P.POSITION,/KEEP_ASPECT_RATIO,true=1
	
	cgimage,bytscl(nir,min=-10,max=80),/keep_aspect,$
	position=!P.POSITION;,MARGIN=0,/MINUS_ONE

	; cgimage,alog(iras),/keep_aspect,$
	; position=!P.POSITION,MARGIN=0,/MINUS_ONE
	
	; print,sp_x[3],sp_y[3],sp_name[3]
	
	imcontour,nir,nirH,/NOERASE, TYPE=0,$
	xmid=sp_x[3],ymid=sp_y[3],$
	color="black",xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" "
	
	
	
	; imcontour,iras,irasH,/NOERASE, TYPE=0,$
	; XSTYLE=1, YSTYLE=1,color=255,xticklen=0.00001,yticklen=0.00001,/nodata,$
	; position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio


	avmapLevels=[2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1]
	
	print,max(avmap),min(avmap)
	; avmaplevelnumber=10
	; avmapMin=max(avmap)*0.500
	; avmapStep=max(avmap)*0.05
	; avmapLevels=findgen(avmaplevelnumber)*avmapStep + avmapMin
	cgcontour,avmap,$
	levels=avmapLevels,Color='black',c_linestyle=0,c_thick=[1,1,1,1,1,1,1,1],$
	c_labels=[0,1,0,1,0,1,0,1],$
	position=!P.POSITION,/OVERPLOT
	; 
	; 
	; 
	; comaplevelnumber=8
	; comapMin=max(comap)*0.200
	; comapStep=max(comap)*0.100
	; comapLevels=findgen(comaplevelnumber)*comapStep + comapMin
	; 
	; cgcontour,comap,$
	; levels=comapLevels,Color='grey',c_linestyle=0,c_thick=[1,1,1,1],$
	; c_labels=[0,0,0,0,0,0,0,0,0,0,0,0],$
	; position=!P.POSITION,/OVERPLOT

cgLoadct,13, CLIP=[185,255],/Silent

cgcolorbar,MINRANGE=-1.85,MAXRANGE=2.74,divisions=5,format='(F0.2)',$
charsize=0.8,Position=[0.6,0.95,0.96,0.97],/TOP

cgplots,sp_x,sp_y,psym=16,symsize=1.1,color=BytScl(sp_alpha),thick=2
cgplots,sp_x,sp_y,psym=9,symsize=1.2,color='black',thick=2

; 
; tvcircle,18,sp_x,sp_y,color=BytScl(sp_alpha),thick=2,/DATA,/FILL
; xyouts,sp_x,sp_y,sp_name,charsize=0.7,alignment=-0.4,color="white"

; for p=0, last_sp-1 do begin
; 	if (((sp_x[p] ge 0) and (sp_x[p] le naxis1)) and ((sp_y[p] ge 0) and (sp_y[p] le naxis2))) then begin
; 		if (sp_class[p] lt 2) then begin
; 			tvcircle,0.1,sp_x[p],sp_y[p],color="red",thick=6,/DATA
; 			xyouts,sp_x[p],sp_y[p],sp_name[p],charsize=0.5,alignment=0.3,color="white"
; 		endif else begin
; 			tvcircle,0.1,sp_x[p],sp_y[p],color="orange",thick=3,/DATA
; 			xyouts,sp_x[p],sp_y[p],sp_name[p],charsize=0.5,alignment=0.3,color="white"
; 			; plots,[sp_x[p],sp_y[p]],psym=8,color=cyan,symsize=2.5,thick=3,/DATA
; 		endelse
; 	endif
; endfor

cgLoadct,0, /Reverse
	
	device,/close
	set_plot,'x'
	!p.font=-1 
	; !P.MULTI=0;
end