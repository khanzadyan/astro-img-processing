pro h2only
	; This program intended to make a figure of H2 line observations
	nb2122=readfits('./NIR/nb2122wcs.fits',nb2122H)
	nb2144=readfits('./NIR/nb2144wcs.fits',nb2144H)
	
	hastrom,nb2144,nb2144H,nb2122H,MISSING=0
	
	; subtraction of continuum from line data
	;
	sub=(nb2122 - nb2144)
	subtr=leefilt(sub,1,1,/EXACT,/DOUBLE)
	subtr(where(subtr le 0)) = -!VALUES.F_NAN
	subtr(where(subtr ge 219)) = !VALUES.F_NAN
	
	naxis1 = sxpar(nb2122H, 'NAXIS1')
	naxis2 = sxpar(nb2122H, 'NAXIS2')
	
	;writefits,'contsub.fits',subtr,nb2122H
	
	FMT = 'A,I,I,F,I,I,F'; defining the format for the ascii file
	extast,nb2122H,h2img_astrometry ; extracting the astrometry info
	
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
	
	READCOL,'catalogs/mhocat.txt',F=FMT,$
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

	
; Part for calculating the outflow sizes

dis=1170 ; parsecs

; part for MHO 735 outflow
print,'=================== part for MHO 735 flow =============================='
print,''
fl1=0.45*SQRT((ABS(h2obj_x[7] - star_x[2]))^2 + (ABS(h2obj_y[7] - star_y[2]))^2)
print,'MHO 735 outflow from IRS3 is',fl1,' arcseconds long'
print,'the real size would then be ',(fl1*dis)/206265,' pc long'
print,''
print,''
;==========================
print,'=================== part for MHO 737 flow =============================='
print,''
fl2=0.45*SQRT((ABS(h2obj_x[9] - star_x[5]))^2 + (ABS(h2obj_y[9] - star_y[5]))^2)
print,'MHO 737 outflow from IRS6 is',fl2,' arcseconds long'
print,'the real size would then be ',(fl2*dis)/206265,' pc long'
print,''
print,''
;==========================
print,'=================== part for 738 and 741 flow =============================='
print,''
fl3p1=SQRT((ABS(h2obj_x[10] - star_x[3]))^2 + (ABS(h2obj_y[10] - star_y[3]))^2)
print, 'part with ',h2obj_name[10],' and ',star_name[3],' = ',fl3p1*0.45
fl3p2=SQRT((ABS(h2obj_x[16] - star_x[3]))^2 + (ABS(h2obj_y[16] - star_y[3]))^2)
print, 'part with ',h2obj_name[16],' and ',star_name[3],' = ',fl3p2*0.45
fl3=0.45*(fl3p1+fl3p2)
print,'together 738 to 741 flow is',fl3,' arcseconds long'
print,'the real size would then be ',(fl3*dis)/206265,' pc long'
print,''
print,''
;==========================
print,'=================== part for 739 and 744 flow =============================='
print,''
fl4p1=SQRT((ABS(h2obj_x[13] - star_x[3]))^2 + (ABS(h2obj_y[13] - star_y[3]))^2)
print, 'part with ',h2obj_name[13],' and ',star_name[3],' = ',fl4p1*0.45
fl4p2=SQRT((ABS(h2obj_x[22] - star_x[3]))^2 + (ABS(h2obj_y[22] - star_y[3]))^2)
print, 'part with ',h2obj_name[22],' and ',star_name[3],' = ',fl4p2*0.45
fl4=0.45*(fl4p1+fl4p2)
print,'together 739 to 744 flow is',fl4,' arcseconds long'
print,'the real size would then be ',(fl4*dis)/206265,' pc long'
print,''
print,''
;=========================
print,'========================== part for 745 flow =============================='
print,''
fl5=0.45*SQRT((ABS(h2obj_x[24] - h2obj_x[26]))^2 + (ABS(h2obj_y[24] - h2obj_y[26]))^2)
print,' MHO 745 outflow is ',fl5,' arcsecond long'
print,'the real size would then be ',(fl5*dis)/206265,' pc long'
print,''
print,''
;=========================
print,'========================== part for 734 and 740 flow ======================'
print,''
fl6=0.45*SQRT((ABS(h2obj_x[0] - h2obj_x[14]))^2 + (ABS(h2obj_y[0] - h2obj_y[14]))^2)
print,' MHO 734 to 740 outflow is ',fl6,' arcsecond long'
print,'the real size would then be ',(fl6*dis)/206265,' pc long'
print,''
print,''
;=========================
print,'=================== part for 742 and 743 flow =============================='
print,''
fl7p1=SQRT((ABS(h2obj_x[18] - star_x[6]))^2 + (ABS(h2obj_y[18] - star_y[6]))^2)
print, 'part with ',h2obj_name[18],' and ',star_name[6],' = ',fl7p1*0.45
fl7p2=SQRT((ABS(h2obj_x[21] - star_x[6]))^2 + (ABS(h2obj_y[21] - star_y[6]))^2)
print, 'part with ',h2obj_name[21],' and ',star_name[6],' = ',fl7p2*0.45
fl7=0.45*(fl7p1+fl7p2)
print,'together 742 to 743 flow is',fl7,' arcseconds long'
print,'the real size would then be ',(fl7*dis)/206265,' pc long'
print,''
print,''


;   
; ======================== Ploting parameters ===================
	
	!X.MARGIN=0
	!Y.MARGIN=0
	!X.THICK=3
	!Y.THICK=3
	!P.CHARSIZE=-1

	loadct,0
	!P.MULTI=[0,1,1]
	!P.POSITION=[0.15,0.12,0.99,0.96]
	!p.font=0 ; nice post script fonts
	set_plot,'ps',/copy, /interpolate
	device,filename='./results/h2only.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
		xsize=15.912,ysize=17,scale_factor=1.0
;   ; ===============================================================
;   
   ;!P.POSITION=[0.10,0.10,0.98,0.98]
   reverse_colors
   tvimage,bytscl(subtr,min=0,max=46,/NAN),$
   	/keep_aspect_ratio,position=!P.POSITION
;   
   reverse_colors
   imcontour,subtr,nb2122H,/NOERASE, TYPE=1,$
   xdelta=2,ydelta=1,$
   ;xmid=star_x(1),ymid=star_y(1),$
   charsize=1,$ xthick=3,ythick=3,$
   XSTYLE=1, YSTYLE=1,color=1,xticklen=0.015,yticklen=0.015,/nodata,$
   position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
;
	@colors
	for s=0, last_star-1 do begin
		if (((star_x[s] ge 0) and (star_x[s] le naxis1)) and ((star_y[s] ge 0) and (star_y[s] le naxis2))) then begin
			tvcross,14,star_x[s],star_y[s],color=red,thick=6,/DATA
		endif
	endfor
;
;   ; H2 knots
	for h=0, last_h2obj-1 do begin
		tvbox,12,h2obj_x[h],h2obj_y[h],color=dkgreen,thick=2,/DATA
	endfor
;   
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   loadct,0
   device,/close
   set_plot,'x'
   !p.font=-1 ; reverting to the IDL default font
   !P.MULTI=0
end