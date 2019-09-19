;################################################################
; This program is to iteratively read in the mcrx stellar
; velocities and fit them first with a varying PA,
; then, holding the PA fixed, fit them again and determine
; both asymmetry parameters (v and sig).
;
; The outputs for this program are a file 'myfile_myr_view.txt'
; that shows the circular velocity of the best fit model, and
; a file text_out_myr_view.txt, which holds the kinematic best
; fitting PA and the asymmetry parameters.
;
; Becky Nevin
;###############################################################

PRO run_kinemetry_iteratively_8318-12704

plateid = '8318-12704'
file = 'kinemetry_input_8318-12704.txt'
file = STRCOMPRESS(file, /REMOVE_ALL)
result = FILE_TEST(file)
IF result EQ 0 THEN CONTINUE
rdfloat, file, xcen, ycen, x_old, y_old, SKIPLIN=1
xcen0=xcen[0]-x_old[0]
ycen0=ycen[0]-y_old[0]
rdfloat, file, num, xbin, ybin, velbin, er_velbin, sigbin, er_sigbin, SKIPLIN=2
range = where(ABS(er_velbin/velbin) lt 1 and ABS(er_sigbin/sigbin) lt 1)
velbin=velbin[range]
er_velbin=er_velbin[range]
xbin=xbin[range]
ybin=ybin[range]
sigbin=sigbin[range]
er_sigbin=er_sigbin[range]
KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, X0=xcen, Y0=ycen, ntrm=10, ERROR=er_velbin, name='simulation',er_cf=er_cf, er_pa=er_pa, er_q=er_q, /verbose, /ALL, /plot, /FIXCEN
k0 = cf[*,0]
k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
k5 = SQRT(cf[*,5]^2 + cf[*,6]^2)
k51 = k5/k1
erk1 = (SQRT( (cf[*,1]*er_cf[*,1])^2 + (cf[*,2]*er_cf[*,2])^2 ))/k1
erk5 = (SQRT( (cf[*,5]*er_cf[*,5])^2 + (cf[*,6]*er_cf[*,6])^2 ))/k5
erk51 = ( SQRT( ((k5/k1) * erk1)^2 + erk5^2  ) )/k1
PRINT, 'median PA before', MEDIAN(pa)-90
index = WHERE(pa GT 180)
pa[index] = pa[index]-180
KIN_PA = MEDIAN(pa)-90
PRINT, 'median PA', KIN_PA
KIN_PA_e = MEDIAN(er_pa/pa)
PA_before = MEDIAN(pa)-90
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
	PRINT, 'Triangulate thing: ', Error_status
	PRINT, !ERROR_STATE.MSG
	CATCH, /CANCEL
	CONTINUE
ENDIF
KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, X0=xcen,Y0=ycen, ntrm=10, ERROR=er_velbin, name='simulation',RANGEPA=[MEDIAN(pa)-90-1,MEDIAN(pa)-90+1],er_cf=er_cf, er_pa=er_q=er_q,/ALL, /FIXCEN, /plot
PAQ1 = [PA_before, MEDIAN(q)]
PAQ = Reform(Rebin(PAQ1, 2, 100), 200)
PRINT, 'about to do only one term'
KINEMETRY, xbin, ybin, velbin, rad, pa, q, cf, X0=xcen,Y0=ycen, ntrm=10, ERROR=er_velbin,  PAQ=PAQ,er_cf=er_cf, er_pa=er_pa, er_q=er_q,velcirc=velCirc,  velkin=velKin, /ALL, /BMODEL, /FIXCEN, /plot
print, cf
k0 = cf[*,0]
k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
k2 = SQRT(cf[*,3]^2 + cf[*,4]^2)
k3 = SQRT(cf[*,5]^2 + cf[*,6]^2)
k4 = SQRT(cf[*,7]^2 + cf[*,8]^2)
k5 = SQRT(cf[*,9]^2 + cf[*,10]^2)
PRINT, 'k0', k0
elemental_v_asym=dblarr(N_ELEMENTS(pa))
FOR I=0, N_ELEMENTS(pa) - 1 DO BEGIN
	k2 = SQRT(cf[i,3]^2 + cf[i,4]^2)
	k3 = SQRT(cf[i,5]^2 + cf[i,6]^2)
	k4 = SQRT(cf[i,7]^2 + cf[i,8]^2)
	k5 = SQRT(cf[i,9]^2 + cf[i,10]^2)
	sum = (k2 + k3 + k4 + k5)/4
	B_1_v = ABS(cf[i,2])
	elemental_v_asym[i]=sum/B_1_v
ENDFOR
alt_v_asym = MEDIAN(elemental_v_asym)
KINEMETRY, xbin, ybin, sigbin, rad, pa, q, cf, X0=xcen, Y0=ycen,ntrm=10, scale=0.5, ERROR=er_sigbin, name='sigma map', PAQ=PAQ, er_cf=er_cf, er_pa=er_pa, er_q=er_q, /verbose, /ALL, /plot, /FIXCEN
k0 = cf[*,0]
k1 = SQRT(cf[*,1]^2 + cf[*,2]^2)
k2 = SQRT(cf[*,3]^2 + cf[*,4]^2)
k3 = SQRT(cf[*,5]^2 + cf[*,6]^2)
k4 = SQRT(cf[*,7]^2 + cf[*,8]^2)
k5 = SQRT(cf[*,9]^2 + cf[*,10]^2)
elemental_s_asym=dblarr(N_ELEMENTS(pa))
FOR I=0, N_ELEMENTS(pa) - 1 DO BEGIN
	k1 = SQRT(cf[i,1]^2 + cf[i,2]^2)
	k2 = SQRT(cf[i,3]^2 + cf[i,4]^2)
	k3 = SQRT(cf[i,5]^2 + cf[i,6]^2)
	k4 = SQRT(cf[i,7]^2 + cf[i,8]^2)
	k5 = SQRT(cf[i,9]^2 + cf[i,10]^2)
	sum = (k1 + k2 + k3 + k4 + k5)/5
	A_0_v = ABS(cf[i,0])
	elemental_s_asym[i]=sum/A_0_v
ENDFOR
alt_s_asym = MEDIAN(elemental_s_asym)
file_out='myfile_v_8318-12704.txt'
file_out = STRCOMPRESS(file_out, /REMOVE_ALL)
writecol, file_out, xbin, ybin, velbin, er_velbin, velCirc, velKin, fmt='(f,f,f,f,x,f,x,f)'
file_out_2 = 'text_out_8318-12704.txt'
file_out_2 = STRCOMPRESS(file_out_2, /REMOVE_ALL)
K_asym=SQRT(alt_s_asym^2+alt_v_asym^2)
writecol, file_out_2, KIN_PA, KIN_PA_e, alt_v_asym, alt_s_asym, K_asym
PRINT, 'v_asym', alt_v_asym, 's_asym', alt_s_asym, 'K_TOT', K_asym
END
