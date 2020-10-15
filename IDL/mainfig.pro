pro mainfig,pic

; This program is for creating sub-figures for the main image in the article 
; dedicated to GM 2-4 region. The main idea is to use switches to create an
; image for different panels.

;
; This is a key to produce a specific figure with this program
; the options are 0=h2, 1=s2 or 2=rgb spitzer

if not keyword_set( pic ) then pic = 0
	
; Reading the input files here ============================================
	
	nb2122=readfits('./NIR/nb2122wcs.fits',nb2122H)
	S2=readfits('./Optical/Byurakan/SII_wcs_rot.fits',S2H)
	
	RChan=readfits('./Infrared/ch4-pbcd-al.fits',h1)
	GChan=readfits('./Infrared/ch2-pbcd-al.fits',h2)
	BChan=readfits('./Infrared/ch1-pbcd-al.fits',h3)
	
; WCS for each input is remapped according to the first input image
	
	hastrom,S2,S2H,nb2122H,MISSING=0
	hastrom,RChan,h1,nb2122H,MISSING=0
	hastrom,GChan,h2,nb2122H,MISSING=0
	hastrom,BChan,h3,nb2122H,MISSING=0
; ====================================== Finish reading files =============

; Scaling of the field ====================================================

	naxis1 = sxpar(nb2122H, 'NAXIS1')
	naxis2 = sxpar(nb2122H, 'NAXIS2')
	
	print,naxis1,naxis2
	
	obj_x=(naxis1-1)/2
	obj_y=((naxis2-1)/2)+36
	
	Xdim=naxis1-1
	Ydim=naxis2-1
	print,Xdim,Ydim
	obj_x=Xdim/2
	obj_y=Ydim/2
	
	Start_x=obj_x - Xdim/2
	End_x=obj_x + Xdim/2
	Start_y=obj_y - Ydim/2
	End_y=obj_y + Ydim/2
; =======================================Finish scaling ===================

; Extracting the sub-fields for each filter ===============================

; H2 ---

	h2img=nb2122(Start_x:End_x,Start_y:End_y)
	hextract,nb2122,nb2122H,h2img,h2imgH,Start_x,End_x,Start_y,End_y

; S2 ---
	s2img=S2(Start_x:End_x,Start_y:End_y)
	hextract,S2,S2H,s2img,s2imgH,Start_x,End_x,Start_y,End_y

; Red channel ---
	Red=RChan(Start_x:End_x,Start_y:End_y)
	hextract,RChan,h1,Red,RedH,Start_x,End_x,Start_y,End_y
	Redflt=leefilt(Red,5,5,/DOUBLE)
	Red(where(Red le 350)) = -!VALUES.F_NAN
	Red(where(Red ge 500)) = !VALUES.F_NAN
	Red=alog10(Red)
	help,Red

; Green channel ---
	Green=GChan(Start_x:End_x,Start_y:End_y)
	hextract,GChan,h2,Green,GreenH, Start_x,End_x,Start_y,End_y
	Greenflt=leefilt(Green,5,5,/DOUBLE)
	Green(where(Green le 13)) = -!VALUES.F_NAN
	Green(where(Green ge 77)) = !VALUES.F_NAN
	Green=alog10(Green)

; Blue channel ---
	Blue=BChan(Start_x:End_x,Start_y:End_y)
	hextract,BChan,h3,Blue,BlueH, Start_x,End_x,Start_y,End_y
	Blueflt=leefilt(Blue,5,5,/DOUBLE)
	Blue(where(Blue le 12)) = -!VALUES.F_NAN
	Blue(where(Blue ge 100)) = !VALUES.F_NAN
	Blue=alog10(Blue)
; ================================== End of Sub-field extraction =========

; Reading the catalougs ==================================================

	FMT = 'A,I,I,F,I,I,F'; defining the format for the ascii file
	extast,h2imgH,h2img_astrometry ; extracting the astrometry info

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

	ad2xy, star_ra, star_dec, h2img_astrometry, star_x, star_y

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
	
	ad2xy, h2obj_ra, h2obj_dec, h2img_astrometry, h2obj_x, h2obj_y

; ============================= End of Catalouge reading ===================

; Loop chooser =============================================================
	loadct,0
	!P.MULTI=[0,1,1]
	!P.POSITION=[0.15,0.12,0.99,0.96]
	!p.font=0 ; nice post script fonts
	set_plot,'ps',/copy, /interpolate

	if pic EQ 0 then begin
		print, 'Will Produce the H2 Figure'
		goto,h2loop
	endif 	
	if pic EQ 1 then begin
		print, 'Will Produce the S2 Figure'
		goto,s2loop
	endif
	if pic EQ 2 then begin
		print, 'Will Produce the Spitzer RGB Figure'
		goto,rgbloop
	endif

h2loop: ; =============== H2LOOP =================================

	device,filename='results/h2mainfig.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
		xsize=15.912,ysize=17,scale_factor=1.0

	reverse_colors
	tvimage,bytscl(h2img,min=-10,max=46),$
		/keep_aspect_ratio,position=!P.POSITION

	reverse_colors
	@colors
	imcontour,h2img,h2imgH,/NOERASE, TYPE=1,$
		xdelta=2,ydelta=1,xthick=3,ythick=3,$
		XSTYLE=1, YSTYLE=1,color=black,xticklen=0.015,yticklen=0.015,/nodata,$
		position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	; h2 objects
	for h=0, last_h2obj-1 do begin
	   tvbox,12,h2obj_x[h],h2obj_y[h],color=dkgreen,thick=4,/DATA
	endfor

	;Potential Sources for Outflows
for s=0, last_star-1 do begin
	if ((star_x[s] ge 0) and (star_x[s] le Xdim) and (star_y[s] ge 0) and (star_y[s] le Ydim+1)) then begin
		tvcross,12,star_x[s],star_y[s],color=red,thick=6,/DATA
	endif
endfor

	loadct,0
    device,/close
	goto,loopend
; ========== End of H2 Loop ======================================

s2loop: ; =============== s2 Loop ================================

	device,filename='results/s2mainfig.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
		xsize=15.912,ysize=17,scale_factor=1.0

	reverse_colors
	tvimage,bytscl(s2img,min=-5,max=60),$
		/keep_aspect_ratio,position=!P.POSITION

	reverse_colors
	@colors
	imcontour,s2img,s2imgH,/NOERASE, TYPE=1,$
		xdelta=2,ydelta=1,xthick=3,ythick=3,$
		XSTYLE=1, YSTYLE=1,color=black,xticklen=0.015,yticklen=0.015,/nodata,$
		position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio

	;Potential Sources for Outflows
	for s=0, last_star-1 do begin
		if ((star_x[s] ge 0) and (star_x[s] le Xdim) and (star_y[s] ge 0) and (star_y[s] le Ydim+1)) then begin
			tvcross,12,star_x[s],star_y[s],color=red,thick=6,/DATA
		endif
	endfor

	loadct,0
    device,/close
	goto,loopend
; ========== End of S2 Loop ======================================

rgbloop: ; ============= Spitzer RGB figure Loop =================

	device,filename='results/rgbmainfig.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
		xsize=15.912,ysize=17,scale_factor=1.0

	axes=lonarr(3) 
	final_image=lonarr(3,Xdim,Ydim+1)
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
	; h2 objects
	for h=0, last_h2obj-1 do begin
		tvbox,12,h2obj_x[h],h2obj_y[h],color=yellow,thick=4,/DATA
	endfor

	;Potential Sources for Outflows
	for s=0, last_star-1 do begin
		if ((star_x[s] ge 0) and (star_x[s] le Xdim) and (star_y[s] ge 0) and (star_y[s] le Ydim+1)) then begin
			tvcross,12,star_x[s],star_y[s],color=red,thick=6,/DATA
		endif
	endfor

	loadct,0
    device,/close
	goto,loopend
; ========== End of RGB Loop ======================================

loopend:
   set_plot,'x'
   !p.font=-1 ; reverting to the IDL default font
   !P.MULTI=0;
end