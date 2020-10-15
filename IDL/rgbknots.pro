pro rgbknots

	; This program is for creating an RGB figure for GM 2-4 region
	
	; Reading the input files here
	
	nb2122=readfits('./NIR/nb2122wcs.fits',nb2122H)
	nb2144=readfits('./NIR/nb2144wcs.fits',nb2144H)
	
	RChan=readfits('./Infrared/ch4-pbcd-al.fits',h1)
	GChan=readfits('./Infrared/ch2-pbcd-al.fits',h2)
	BChan=readfits('./Infrared/ch1-pbcd-al.fits',h3)
	
	; WCS for each input is remapped according to the first input image
	
	hastrom,nb2144,nb2144H,nb2122H,MISSING=0
	hastrom,RChan,h1,nb2122H,MISSING=0
	hastrom,GChan,h2,nb2122H,MISSING=0
	hastrom,BChan,h3,nb2122H,MISSING=0

	naxis1 = sxpar(nb2122H, 'NAXIS1')
	naxis2 = sxpar(nb2122H, 'NAXIS2')
	
	Xdim=naxis1-1
	Ydim=naxis2-1
	
	obj_x=Xdim/2
	obj_y=Ydim/2
	
	Start_x=obj_x - Xdim/2
	End_x=obj_x + Xdim/2
	Start_y=obj_y - Ydim/2
	End_y=obj_y + Ydim/2

	; These are the coordinates for the object in interest
	;
	; 05:40:24.3 23:50:55.2 our object
	;obj_ra=15*ten(05,40,24.3)
	;obj_dec=ten(23,50,55.2) 
	
	; 05:40:27.568 +23:52:17.62 HH941A
	;obj_ra=15*ten(05,40,27.6)
	;obj_dec=ten(23,52,17.6)
	;
	; Ploting parameters ----
	loadct,0	
	!P.MULTI=[0,1,1]
	!p.font=0
	set_plot,'ps',/copy, /interpolate
	device,filename='spitzer3col-knots.eps',/encapsulated,/color,bits=8, $
		/helvetica,/isolatin1,$
		xoffset=0,yoffset=0,$
		xsize=11*1.2,ysize=11*1.2,scale_factor=1.0
	
	
	PicPos=[0.10,0.10,0.98,0.98]
	; -----------------------
		
	; ================= PART FOR THE WIDE FIELD IMAGE =====================
	; Here we determine the size of a box centered on the object
	
		
	;adxy,h1,obj_ra,obj_dec,obj_x,obj_y
	
	;HalfLength=100
	;Start_x=obj_x - HalfLength
	;End_x=obj_x + HalfLength
	;Start_y=obj_y - HalfLength
	;End_y=obj_y + HalfLength
		
	Red=RChan(Start_x:End_x,Start_y:End_y)
	hextract,RChan,h1,Red,RedH,Start_x,End_x,Start_y,End_y
	print,min(Red),max(Red)
	Redflt=leefilt(Red,5,5,/DOUBLE)
	Red(where(Red le 350)) = -!VALUES.F_NAN
	Red(where(Red ge 500)) = !VALUES.F_NAN
	print,min(Red),max(Red)
	;maxred=400
	;minred=360
	;maxred=600
	;minred=min(Red)
	Red=alog10(Red)
	
	
	Green=GChan(Start_x:End_x,Start_y:End_y)
	hextract,GChan,h2,Green,GreenH, Start_x,End_x,Start_y,End_y
	print,min(Green),max(Green)
	Greenflt=leefilt(Green,5,5,/DOUBLE)
	Green(where(Green le 13)) = -!VALUES.F_NAN
	Green(where(Green ge 77)) = !VALUES.F_NAN
	print,min(Green),max(Green)
	;maxgreen=77
	;mingreen=14
	;maxgreen=200
	;mingreen=min(Green)
	Green=alog10(Green)
	
	
	Blue=BChan(Start_x:End_x,Start_y:End_y)
	hextract,BChan,h3,Blue,BlueH, Start_x,End_x,Start_y,End_y
	print,min(Blue),max(Blue)
	Blueflt=leefilt(Blue,5,5,/DOUBLE)
	Blue(where(Blue le 12)) = -!VALUES.F_NAN
	Blue(where(Blue ge 100)) = !VALUES.F_NAN
	print,min(Blue),max(Blue)
	;maxblue=51
	;minblue=15
	;maxblue=200
	;minblue=min(Blue)
	Blue=alog10(Blue)
	
	naxis1 = sxpar(RedH, 'NAXIS1')
	naxis2 = sxpar(RedH, 'NAXIS2')
	
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
	
	
	
	
	axes=lonarr(3) 
	final_image=lonarr(3,naxis1,naxis2)
	; construction of rgb image from the suplied min and max
	final_Image[0,*,*] = bytscl(Red,min=minred,max=maxred)
	final_Image[1,*,*] = bytscl(Green,min=mingreen,max=maxgreen)
	final_Image[2,*,*] = bytscl(Blue,min=minblue,max=maxblue)
	
	!P.POSITION=PicPos
	
	tvscale,final_image,position=Pos,/KEEP_ASPECT_RATIO,true=1
	
	imcontour,Red,RedH,/NOERASE, TYPE=1,charsize=1.,$
	xmid=star_x(1),ymid=star_y(1),$
	XSTYLE=1, YSTYLE=1,color=1,xticklen=-0.015,yticklen=-0.015,/nodata,$
	position=Pos, SUBTITLE=" " 
	
	; ===== > Overploting the objects
	
	@colors ; reading custom colours
	
	
	; H2 knots
		
	ad2xy, h2obj_ra, h2obj_dec, h2img_astrometry, h2obj_x, h2obj_y
	
	for h=0, last_h2obj-1 do begin
		tvbox,12,h2obj_x[h],h2obj_y[h],color=yellow,thick=2,/DATA
	endfor
	
	k=0
	for k=0, last_star-1 do begin
			tvcross,10,star_x[k],star_y[k],color=red,thick=6,/DATA
	endfor
	
	device,/close
spawn,'/usr/local/bin/convert -density 300 -quality 100 spitzer3col-knots.eps spitzer3col-knots.png'
spawn,'open spitzer3col-knots.png'
	set_plot,'x'
	!p.font=-1 ; reverting to the IDL default font
	!P.MULTI=0

end
