pro fieldall,pic

	; This program is for creating a figures for the IRS objects
	; This program was used to create separate sub-images used on the first image of the 
    ; article on GM2-4 region.
	
if not keyword_set( pic ) then pic = 0

; Reading the input files here ============================================	
	refimg=readfits('./Optical/dss2-red.fits',refimgH)
	;
	avmapimg=readfits('./Optical/78023543-fk5-rot.fits',avmapimgH)
	;
	Jimg=readfits('./NIR/2MASS-J.fits',JimgH)
	Himg=readfits('./NIR/2MASS-H.fits',HimgH)
	Kimg=readfits('./NIR/2MASS-K.fits',KimgH)
	;
	irac1img=readfits('./Infrared/ch1-pploc-al.fits',irac1imgH)
	irac2img=readfits('./Infrared/ch2-pploc-al.fits',irac2imgH)		
	irac4img=readfits('./Infrared/ch4-pploc-al.fits',irac4imgH)
	;
	mips1img=readfits('./Infrared/24mu-pbcd-al.fits',mips1imgH)
	mips2img=readfits('./Infrared/70mu-pbcd-al.fits',mips2imgH)
	mips3img=readfits('./Infrared/160mu-pbcd-al.fits',mips3imgH)
	;
	iras4img=readfits('./FIR/iras100mu.fits',iras4imgH)
	;
	scub1img=readfits('./Radio/m03au44short.fits',scub1imgH)
	scub2img=readfits('./Radio/m03au44long.fits',scub2imgH)
	;
	cmimg=readfits('./Radio/m012_c1_new.fits',cmimgH)
; ================================================ Finish reading the files	

; Scaling of the field ====================================================
	
	hastrom,avmapimg,avmapimgH,refimgH,MISSING=0
	;
	hastrom,Jimg,JimgH,refimgH,MISSING=0
	hastrom,Himg,HimgH,refimgH,MISSING=0
	hastrom,Kimg,KimgH,refimgH,MISSING=0
	;
	hastrom,irac1img,irac1imgH,refimgH,MISSING=0
	hastrom,irac2img,irac2imgH,refimgH,MISSING=0
	hastrom,irac4img,irac4imgH,refimgH,MISSING=0
	;
	hastrom,mips1img,mips1imgH,refimgH,MISSING=0
	hastrom,mips2img,mips2imgH,refimgH,MISSING=0
	hastrom,mips3img,mips3imgH,refimgH,MISSING=0
	;
	hastrom,iras4img,iras4imgH,refimgH,MISSING=0
	;
	hastrom,scub1img,scub1imgH,refimgH,MISSING=0
	hastrom,scub2img,scub2imgH,refimgH,MISSING=0
	;
	hastrom,cmimg,cmimgH,refimgH,MISSING=0
	

	; 05 40 24.233 +23 50 54.62
	obj_ra=15*ten(05,40,24.233)
	;obj_ra=15*ten(05,40,20.000)
	obj_dec=ten(23,50,54.62)
	
	; defining the center of the field
	adxy,refimgH,obj_ra,obj_dec,obj_x,obj_y
	HalfLength=300
	Start_x=obj_x - HalfLength
	End_x=obj_x + HalfLength
	Start_y=obj_y - HalfLength
	End_y=obj_y + HalfLength
; ======================================= Finish scaling of the field ========

; Extracting the sub-fields for each filter ===============================	
	
	dss=refimg(Start_x:End_x,Start_y:End_y)
	hextract,refimg,refimgH,dss,dssH,Start_x,End_x,Start_y,End_y
	sky,dss,dss_skyvalue,dss_skysig,/NAN
	dss(where(dss lt dss_skyvalue)) = dss_skyvalue
	
	naxis1 = sxpar(dssH, 'NAXIS1')
	naxis2 = sxpar(dssH, 'NAXIS2')
	
	
; AV MAP ===

	avmap=avmapimg(Start_x:End_x,Start_y:End_y)
	hextract,avmapimg,avmapimgH,avmap,avmapH,Start_x,End_x,Start_y,End_y
	sky,avmap,avmap_skyvalue,avmap_skysig,/NAN

; J-band ===
	
	massBlue=Jimg(Start_x:End_x,Start_y:End_y)
	hextract,Jimg,JimgH,massBlue,massBlueH, Start_x,End_x,Start_y,End_y
	sky,massBlue,massBlue_skyvalue,massBlue_skysig,/NAN,/SILENT
	massBlue(where(massBlue lt (massBlue_skyvalue))) = (massBlue_skyvalue)
	massBlue=alog10(massBlue)
	
; H-band ===
	
	massGreen=Himg(Start_x:End_x,Start_y:End_y)
	hextract,Himg,HimgH,massGreen,massGreenH, Start_x,End_x,Start_y,End_y
	sky,massGreen,massGreen_skyvalue,massGreen_skysig,/NAN,/SILENT
	massGreen(where(massGreen lt (massGreen_skyvalue))) = (massGreen_skyvalue)
	massGreen=alog10(massGreen)
	
; K-band ===	
	
	massRed=Kimg(Start_x:End_x,Start_y:End_y)
	hextract,Kimg,KimgH,massRed,massRedH, Start_x,End_x,Start_y,End_y
	sky,massRed,massRed_skyvalue,massRed_skysig,/NAN,/SILENT
	massRed(where(massRed lt (massRed_skyvalue))) = (massRed_skyvalue)
	massRed=alog10(massRed)
	
; IRAC ch1 ===

	Blue=irac1img(Start_x:End_x,Start_y:End_y)
	hextract,irac1img,irac1imgH,Blue,BlueH, Start_x,End_x,Start_y,End_y
	sky,Blue,Blue_skyvalue,Blue_skysig,/NAN;,/SILENT
	print,min(Blue),max(Blue)
	minblue=-2
	maxblue=15
	;Blue(where(Blue lt 17)) = !VALUES.F_NAN
	; Blue(where(Blue lt (Blue_skyvalue))) = (Blue_skyvalue)
	;Blue(where(Blue lt 0.22)) = !VALUES.F_NAN
	;Blue=alog10(Blue)

; IRAC ch2 ===

	Green=irac2img(Start_x:End_x,Start_y:End_y)
	hextract,irac2img,irac2imgH,Green,GreenH, Start_x,End_x,Start_y,End_y
	sky,Green,Green_skyvalue,Green_skysig,/NAN;,/SILENT
	print,min(Green),max(Green)
	mingreen=-2
	maxgreen=15
	;Green(where(Green lt 16)) = !VALUES.F_NAN
	; Green(where(Green lt (Green_skyvalue))) = (Green_skyvalue)
	;Green(where(Green lt 0.3)) = !VALUES.F_NAN
	;Green=alog10(Green)

; IRAC ch4 ===

	Red=irac4img(Start_x:End_x,Start_y:End_y)
	hextract,irac4img,irac4imgH,Red,RedH,Start_x,End_x,Start_y,End_y
	sky,Red,Red_skyvalue,Red_skysig,/NAN;,/SILENT
	print,min(Red),max(Red)
	minred=9
	maxred=15
	;Red(where(Red lt 360)) = !VALUES.F_NAN
	; Red(where(Red lt (Red_skyvalue))) = (Red_skyvalue)
	;Red(where(Red lt 10)) = !VALUES.F_NAN
	;Red=alog10(Red)

; MIPS ch1 ===

	mu24=mips1img(Start_x:End_x,Start_y:End_y)
	hextract,mips1img,mips1imgH,mu24,mu24H,Start_x,End_x,Start_y,End_y
	print,'mu24 part ---------------------------------------'
	sky,mu24,mu24_skyvalue,mu24_skysig,/NAN;,/SILENT
	mu24=alog10(mu24)

; MIPS ch2 ===
	
	mu70=mips2img(Start_x:End_x,Start_y:End_y)
	hextract,mips2img,mips2imgH,mu70,mu70H,Start_x,End_x,Start_y,End_y
	sky,mu70,mu70_skyvalue,mu70_skysig,/NAN,/SILENT
	; mu70(where(mu70 lt (mu70_skyvalue))) = (mu70_skyvalue)
	; mu70(where(mu70 gt 4233)) = !VALUES.F_NAN
	; mu70=alog10(mu70)

; MIPS ch3 ===

	mu160=mips3img(Start_x:End_x,Start_y:End_y)
	hextract,mips3img,mips3imgH,mu160,mu160H,Start_x,End_x,Start_y,End_y
	;mu160=alog10(mu160)

; IRAS ch4 ===

	iras100=iras4img(Start_x:End_x,Start_y:End_y)
	hextract,iras4img,iras4imgH,iras100,iras100H,Start_x,End_x,Start_y,End_y
	sky,iras100,iras100_skyvalue,iras100_skysig,/NAN,/SILENT

; SCUBA 450mu ===

	scuba450mu=scub1img(Start_x:End_x,Start_y:End_y)
	hextract,scub1img,scub1imgH,scuba450,scuba450H,Start_x,End_x,Start_y,End_y
	sky,scuba450,scuba450_skyvalue,scuba450_skysig,/NAN,/SILENT

; SCUBA 850mu ===

	scuba850mu=scub2img(Start_x:End_x,Start_y:End_y)
	hextract,scub2img,scub2imgH,scuba850,scuba850H,Start_x,End_x,Start_y,End_y
	sky,scuba850,scuba850_skyvalue,scuba850_skysig,/NAN,/SILENT
	
; 3.6 cm data ===

	cm=cmimg(Start_x:End_x,Start_y:End_y)
	hextract,cmimg,cmimgH,cm,cmH,Start_x,End_x,Start_y,End_y
	cm(where(cm le 0)) = -!VALUES.F_NAN
	sky,cm,cm_skyvalue,cm_skysig,/NAN,/SILENT
	;cm=alog10(cm)

; =========================================== Finish extracting sub-fields ==

; Reading the catalougs ==================================================

	FMT = 'A,I,I,F,I,I,F'; defining the format for the ascii file
	extast,dssH,dss_astrometry ; extracting the astrometry info

	READCOL,'catalogs/irs.txt',F=FMT,$
	star_name,star_hr,star_min,star_sec,star_deg,star_dmin,star_dsec
	
	starcount=size(star_name)
	last_star=starcount(1)
	star_ra=dindgen(last_star)
	star_dec=dindgen(last_star)
	star_x=fltarr(last_star)
	star_y=fltarr(last_star)

	for s=0, last_star-1 do begin
		star_ra[s]=15*ten(star_hr[s],star_min[s],star_sec[s])
		star_dec[s]=ten(star_deg[s],star_dmin[s],star_dsec[s])
	endfor

	ad2xy, star_ra, star_dec, dss_astrometry, star_x, star_y

	READCOL,'catalogs/mho-final.txt',F=FMT,$
	h2obj_name,h2obj_hr,h2obj_min,h2obj_sec,h2obj_deg,h2obj_dmin,h2obj_dsec

	h2objcount=size(h2obj_name)
	last_h2obj=h2objcount(1)
	h2obj_ra=dindgen(last_h2obj)
	h2obj_dec=dindgen(last_h2obj)	
	h2obj_x=fltarr(last_h2obj)
	h2obj_y=fltarr(last_h2obj)
	
	for h=0, last_h2obj-1 do begin
		h2obj_ra[h]=15*ten(h2obj_hr[h],h2obj_min[h],h2obj_sec[h])
		h2obj_dec[h]=ten(h2obj_deg[h],h2obj_dmin[h],h2obj_dsec[h])
	endfor
	
	ad2xy, h2obj_ra, h2obj_dec, dss_astrometry, h2obj_x, h2obj_y

	FMTSpitzer='A,I,I,F,I,I,F,I'
	READCOL,'catalogs/spitzer.txt',F=FMTSpitzer,$
	sp_name,sp_hr,sp_min,sp_sec,sp_deg,sp_dmin,sp_dsec,sp_class
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

	ad2xy, sp_ra, sp_dec, dss_astrometry, sp_x, sp_y
	; print,sp_name,sp_x,sp_y


; ============================= End of Catalouge reading ===================

; Loop chooser =============================================================
	loadct,0
	!P.MULTI=[0,1,1]
	!P.POSITION=[0.15,0.12,0.99,0.96]
	!p.font=0 ; nice post script fonts
	set_plot,'ps',/copy, /interpolate

	if pic EQ 0 then begin
		print, 'Will Produce the DSS2 Figure'
		goto,ibandloop
	endif 	
	if pic EQ 1 then begin
		print, 'Will Produce the 2MASS RGB Figure'
		goto,massloop
	endif
	if pic EQ 2 then begin
		print, 'Will Produce the Spitzer RGB Figure'
		goto,spitzerloop
	endif
	if pic EQ 3 then begin
		print, 'Will Produce the Spitzer 24mu Figure'
		goto,spitzer24loop
	endif
	if pic EQ 4 then begin
		print, 'Will Produce the Spitzer 70mu Figure'
		goto,spitzer70loop
	endif
	if pic EQ 5 then begin
		print, 'Will Produce the Spitzer 160mu Figure'
		goto,spitzer160loop
	endif
	if pic EQ 6 then begin
		print, 'Will Produce the SCUBA 450mu Figure'
		goto,scuba450loop
	endif
	if pic EQ 7 then begin
		print, 'Will Produce the SCUBA 850mu Figure'
		goto,scuba850loop
	endif
	if pic EQ 8 then begin
		print, 'Will Produce the 3.6cm Figure'
		goto,cmloop
	endif

ibandloop: ; =============== DSS =================================

	device,filename='results/regdss.eps',/encapsulated,/color,bits=8, $
	/helvetica,/isolatin1,$
	xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0

	reverse_colors
	tvimage,bytscl(dss,min=4300,max=8554),$
		/keep_aspect_ratio,position=!P.POSITION

	reverse_colors
	@colors
	imcontour,dss,dssH,/NOERASE, TYPE=1,$
		xdelta=2,ydelta=1,xthick=3,ythick=3,$
		XSTYLE=1, YSTYLE=1,color=black,xticklen=0.015,yticklen=0.015,/nodata,$
		position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio


	avlevels=[2.0,2.2,2.4,2.6,2.8,3.0]
	avlevthick=[1,2,3,4,5,6]
	contour,avmap,$
	levels=avlevels,C_Colors=grey,c_linestyle=0,c_thick=avlevthick,$
	c_labels=[1,1,0,1,1,1],$
	position=!P.POSITION,/OVERPLOT

	; SCUBA 450 -------------
	; Scuba450levelnumber=round(((max(scuba450)-min(scuba450))/scuba450_skysig))
	; Scuba450Min=scuba450_skysig*5
	; Scuba450Step=scuba450_skysig*3
	; Scuba450Levels=findgen(Scuba450levelnumber)*Scuba450Step + Scuba450Min
	
;	Scuba450Levels=[2694,3593,4491,5389,6288,7186,8084,8983,9881,10779,11678,12576,13474,14373,15271,16169,17068,17966]
	; Scuba450Levels=[2694,4491,6288,8084,9881,11678,13474,15271,17068,17966]
	; 
	; contour,scuba450,$
	; levels=Scuba450Levels,C_Colors=blue,c_linestyle=0,c_thick=1,$
	; position=!P.POSITION,/OVERPLOT

	;Potential Sources for Outflows
	; for s=0, last_star-1 do begin
	; 	tvcross,12,star_x[s],star_y[s],color=red,thick=6,/DATA
	; endfor

	; IRAS source
	iras_ra=15*ten(05,40,24.4)
	iras_dec=ten(23,50,54.0)
	adxy,dssH,iras_ra,iras_dec,iras_x,iras_y
	tvellipse,30,6,iras_x,iras_y,88+90,$
	color=blue,thick=4,linestyle=0,/DATA
	
	; maser 1
	mas1_ra=15*ten(05,40,24.9)
	mas1_dec=ten(23,50,56.0)
	adxy,dssH,mas1_ra,mas1_dec,mas1_x,mas1_y
	plots,mas1_x,mas1_y,psym=2,symsize=1,thick=2,color=dkgreen,/DATA
	
	; maser 2
	mas2_ra=15*ten(05,40,24.4)
	mas2_dec=ten(23,50,54.0)
	adxy,dssH,mas2_ra,mas2_dec,mas2_x,mas2_y
	plots,mas2_x,mas2_y,psym=2,symsize=1,thick=2,color=dkgreen,/DATA

	loadct,0
    device,/close
	goto,loopend
; ===================== I band ======================================


massloop: ; =============== 2MASS ==================================

	device,filename='results/reg2massRGB.eps',/encapsulated,/color,bits=8, $
	/helvetica,/isolatin1,$
	xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	axes=lonarr(3) 
	final_image=lonarr(3,naxis1,naxis2)
	; construction of rgb image from the suplied min and max
	final_Image[0,*,*] = bytscl(massRed);,min=minred,max=maxred)
	final_Image[1,*,*] = bytscl(massGreen);,min=mingreen,max=maxgreen)
	final_Image[2,*,*] = bytscl(massBlue);,min=minblue,max=maxblue)
	
	tvscale,final_image,position=!P.POSITION,/KEEP_ASPECT_RATIO,true=1
	
	imcontour,massRed,massRedH,/NOERASE, TYPE=1,charsize=1.,$
	xdelta=2,ydelta=1,xthick=3,ythick=3,$
	;xmid=star_x(1),ymid=star_y(1),$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE="
	
	imcontour,massRed,massRedH,/NOERASE, TYPE=1,charsize=1.,$
	xdelta=2,ydelta=1,xthick=3,ythick=3,$
	;xmid=star_x(1),ymid=star_y(1),$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" "
	
	@colors
	
	; IRAS 100 micron ----------------
	iras100levelnumber=round(((max(iras100)-min(iras100))/iras100_skysig))
	iras100Min=iras100_skysig
	iras100Step=iras100_skysig
	iras100Levels=findgen(iras100levelnumber)*iras100Step + iras100Min
	print,iras100Levels
	
	iras100Levels=[50,70,90,110]
	iras100levelthick=[2,4,6,8]
	
	contour,iras100,$
	levels=iras100Levels,C_Colors=white,c_linestyle=0,$
	C_thick=iras100levelthick,c_labels=[1,1,1,1],$
	position=!P.POSITION,/OVERPLOT
	
	; @colors
	; for s=0, last_star-1 do begin
	; 	tvcross,15,star_x[s],star_y[s],color=white,thick=6,/DATA
	; endfor
	
	loadct,0
    device,/close
	goto,loopend
; ===================== 2MASS =======================================


spitzerloop: ; =============== SPITZER IRAC =========================
	
	device,filename='results/regIracRGB.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0

	axes=lonarr(3) 
	final_image=lonarr(3,naxis1,naxis2)
	; construction of rgb image from the suplied min and max
	final_Image[0,*,*] = bytscl(Red,min=minred,max=maxred)
	final_Image[1,*,*] = bytscl(Green,min=mingreen,max=maxgreen)
	final_Image[2,*,*] = bytscl(Blue,min=minblue,max=maxblue)

	tvscale,final_image,position=!P.POSITION,/KEEP_ASPECT_RATIO,true=1
	
	imcontour,Red,RedH,/NOERASE, TYPE=1,charsize=1.,$
	xdelta=2,ydelta=1,xthick=3,ythick=3,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio 
	
	imcontour,Red,RedH,/NOERASE, TYPE=1,charsize=1.,$
	xdelta=2,ydelta=1,xthick=3,ythick=3,$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	@colors
; 	; psym=-symcat(34)
for p=0, last_sp-1 do begin
	if (((sp_x[p] ge 0) and (sp_x[p] le naxis1)) and ((sp_y[p] ge 0) and (sp_y[p] le naxis2))) then begin
		if (sp_class[p] lt 2) then begin
			tvcircle,12,sp_x[p],sp_y[p],color=ltblue,thick=6,/DATA
			; xyouts,sp_x[p],sp_y[p],sp_name[p],charsize=0.8,alignment=0.3,color=cyan
		endif else begin
			tvcircle,7,sp_x[p],sp_y[p],color=ltblue,thick=3,/DATA
			; plots,[sp_x[p],sp_y[p]],psym=8,color=cyan,symsize=2.5,thick=3,/DATA
		endelse
	endif
endfor

	;Potential Sources for Outflows
	; for s=0, last_star-1 do begin
	; 	tvcross,12,star_x[s],star_y[s],color=red,thick=6,/DATA
	; endfor
	
	loadct,0
    device,/close
	goto,loopend
; ===================== SPITZER IRAC ================================

spitzer24loop: ; ============= SPITZER 24mu =========================
	
	loadct,15
	device,filename='results/regMips24mu.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	;reverse_colors
	tvimage,bytscl(mu24,min=min(mu24),max=max(mu24)),$
		/keep_aspect_ratio,position=!P.POSITION
	;   
	;reverse_colors
	imcontour,mu24,mu24H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,xthick=3,ythick=3,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	imcontour,mu24,mu24H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,xthick=3,ythick=3,$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
	
	; @colors
	; Scuba450Levels=[2694,4491,6288,8084,9881,11678,13474,15271,17068,17966]
	; 
	; contour,scuba450,$
	; levels=Scuba450Levels,C_Colors=black,c_linestyle=0,c_thick=1,$
	; position=!P.POSITION,/OVERPLOT
	
	; @colors	
	; ;Potential Sources for Outflows
	; for s=0, last_star-1 do begin
	; 	tvcircle,24,star_x[s],star_y[s],color=cyan,thick=3,/DATA
	; endfor

	
	loadct,0
    device,/close
	goto,loopend
; ===================== SPITZER 24mu ================================

spitzer70loop: ; ============= SPITZER 70mu =========================
	
	loadct,15
	device,filename='regMips70mu.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	;reverse_colors
	tvimage,bytscl(mu70,min=25,max=1000),$
		/keep_aspect_ratio,position=!P.POSITION
	;   
	;reverse_colors
	imcontour,mu70,mu70H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	imcontour,mu70,mu70H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
	
	@colors
	Scuba450Levels=[2694,4491,6288,8084,9881,11678,13474,15271,17068,17966]
	
	contour,scuba450,$
	levels=Scuba450Levels,C_Colors=black,c_linestyle=0,c_thick=1,$
	position=!P.POSITION,/OVERPLOT
	

	;Potential Sources for Outflows
	for s=0, last_star-1 do begin
		tvcross,12,star_x[s],star_y[s],color=red,thick=6,/DATA
	endfor
	
	loadct,0
    device,/close
	goto,loopend
; ===================== SPITZER 70mu ================================

spitzer160loop: ; ============ SPITZER 160mu ========================
	
	loadct,15
	device,filename='regMips160mu.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	;reverse_colors
	tvimage,bytscl(mu160,min=1,max=230),$
		/keep_aspect_ratio,position=!P.POSITION
	;   
	;reverse_colors
	imcontour,mu160,mu160H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	imcontour,mu160,mu160H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
	
	@colors
	
	; IRAS 100 micron ----------------
	iras100levelnumber=round(((max(iras100)-min(iras100))/iras100_skysig))
	iras100Min=iras100_skysig*3
	iras100Step=iras100_skysig
	iras100Levels=findgen(iras100levelnumber)*iras100Step + iras100Min
	
	contour,iras100,$
	levels=iras100Levels,C_Colors=red,c_linestyle=1,$
	position=!P.POSITION,/OVERPLOT
	
	; SCUBA 850 -------------
	Scuba850levelnumber=round(((max(scuba850)-min(scuba850))/scuba850_skysig))
	Scuba850Min=scuba850_skysig*5
	Scuba850Step=scuba850_skysig*3
	Scuba850Levels=findgen(Scuba850levelnumber)*Scuba850Step + Scuba850Min
	
;	Scuba450Levels=[2694,3593,4491,5389,6288,7186,8084,8983,9881,10779,11678,12576,13474,14373,15271,16169,17068,17966]
		
	contour,scuba850,$
	levels=Scuba850Levels,C_Colors=white,c_linestyle=0,c_thick=1,$
	position=!P.POSITION,/OVERPLOT
	
	;Potential Sources for Outflows
	for s=0, last_star-1 do begin
		tvcross,12,star_x[s],star_y[s],color=red,thick=6,/DATA
	endfor
	
	loadct,0
    device,/close
	goto,loopend
; ===================== SPITZER 160mu ===============================

scuba450loop: ; ============== SCUBA 450mu ==========================
	
	loadct,15
	device,filename='regScuba450mu.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	;reverse_colors
	tvimage,bytscl(scuba450,min=0.1,max=max(scuba450)),$
		/keep_aspect_ratio,position=!P.POSITION
	;   
	;reverse_colors
	imcontour,scuba450,scuba450H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	imcontour,scuba450,scuba450H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
	
	loadct,0
    device,/close
	goto,loopend
; ===================== SCUBA 450mu ==================================

scuba850loop: ; ============== SCUBA 850mu ==========================
	
	loadct,15
	device,filename='regScuba850mu.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	;reverse_colors
	tvimage,bytscl(scuba850,min=0.1,max=max(scuba850)),$
		/keep_aspect_ratio,position=!P.POSITION
	;   
	;reverse_colors
	imcontour,scuba850,scuba850H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	imcontour,scuba850,scuba850H,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
	
	loadct,0
    device,/close
	goto,loopend
; ===================== SCUBA 850mu ==================================

cmloop: ; ============== 3.6cm ====================================
	
	loadct,15
	device,filename='reg3p6cm.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	;reverse_colors
	tvimage,bytscl(cm,min=0.002,max=0.0111526),$
		/keep_aspect_ratio,position=!P.POSITION
	;   
	;reverse_colors
	imcontour,cm,cmH,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	imcontour,cm,cmH,/NOERASE, TYPE=1,$
	charsize=1,$
	xdelta=2,ydelta=1,$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=0.0001,yticklen=0.0001,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
	
	loadct,0
    device,/close
	goto,loopend
; ===================== 3.6cm ========================================



loopend:
   set_plot,'x'
   !p.font=-1 ; reverting to the IDL default font
   !P.MULTI=0;
end