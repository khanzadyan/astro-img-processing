pro field
	; Purpose of this file is to produce picture of whole field in H2
	; over-layed with the known and discovered objects

	dssimg=readfits('./Optical/dss2-red.fits',dssimgH)
	avmapimg=readfits('./Optical/78023543-fk5-rot.fits',avmapimgH)
	
	hastrom,avmapimg,avmapimgH,dssimgH,MISSING=0
	
	
	; 05 40 24.233 +23 50 54.62
	obj_ra=15*ten(05,40,24.233)
	obj_dec=ten(23,50,54.62)
	
	; defining the center of the field
	adxy,dssimgH,obj_ra,obj_dec,obj_x,obj_y
	HalfLength=300
	Start_x=obj_x - HalfLength
	End_x=obj_x + HalfLength
	Start_y=obj_y - HalfLength
	End_y=obj_y + HalfLength
	
	
	dss=dssimg(Start_x:End_x,Start_y:End_y)
	hextract,dssimg,dssimgH,dss,dssH,Start_x,End_x,Start_y,End_y
	sky,dss,dss_skyvalue,dss_skysig,/NAN
	dss(where(dss lt dss_skyvalue)) = dss_skyvalue
	
	
	avmap=avmapimg(Start_x:End_x,Start_y:End_y)
	hextract,avmapimg,avmapimgH,avmap,avmapH,Start_x,End_x,Start_y,End_y
	sky,avmap,avmap_skyvalue,avmap_skysig,/NAN
	
	sizedss=size(dss)
	sizex=sizedss(1)
	sizey=sizedss(2)
		
	loadct,0
	!P.MULTI=[0,1,1]
	!P.POSITION=[0.12,0.12,0.98,0.98]
	!p.font=0 ; nice post script fonts
	set_plot,'ps',/copy, /interpolate
	
	device,filename='whole.eps',/encapsulated,/color,bits=8, $
	/helvetica,/isolatin1,$
	xoffset=0,yoffset=0,$
	xsize=15,ysize=14.97,scale_factor=1.0
	
	reverse_colors
	tvimage,bytscl(dss,min=4300,max=8554),$
	/keep_aspect_ratio,position=!P.POSITION
		
	imcontour,dss,dssH,/NOERASE, TYPE=1,$
	;xdelta=2,ydelta=1,$
	charsize=1,$
	XSTYLE=1, YSTYLE=1,color=255,xticklen=0.015,yticklen=0.015,/nodata,$
	position=!P.POSITION, SUBTITLE=" ",/keep_aspect_ratio
	
	@colors

	avlevels=[2.0,2.2,2.4,2.6,2.8,3.0]
	avlevthick=[1,2,3,4,5,6]
	contour,avmap,$
	levels=avlevels,C_Colors=grey,c_linestyle=0,c_thick=avlevthick,$
	c_labels=[1,1,0,1,1,1],$
	position=!P.POSITION,/OVERPLOT
	
	; IRAS source
	iras_ra=15*ten(05,40,24.4)
	iras_dec=ten(23,50,54.0)
	adxy,dssH,iras_ra,iras_dec,iras_x,iras_y
	tvellipse,30,6,iras_x,iras_y,88+90,$
	color=blue,thick=4,linestyle=0,/DATA
	
	; GGD 4
	; ggd4_ra=15*ten(05,40,25.0)
	; ggd4_dec=ten(23,50,54.0)
	; adxy,dssH,ggd4_ra,ggd4_dec,ggd4_x,ggd4_y
	; tvellipse,18,18,ggd4_x,ggd4_y,98+90,$
	; color=red,thick=4,linestyle=2,/DATA
	
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
	
	device,/close
	set_plot,'x'
	!p.font=-1 
	!P.MULTI=0;
end