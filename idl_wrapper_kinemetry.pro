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

PRO idl_wrapper_kinemetry, list_plateifus

DEVICE, DECOMPOSED=0
FOR j=0, N_ELEMENTS(list_plateifus)-1 DO BEGIN
;	print, STRJOIN(STRSPLIT(STRING(list_plateifus[j]),'-', /EXTRACT), '_')
	RUN_KINEMETRY_ON_MULTIPLE_GALAXIES, list_plateifus[j]
ENDFOR
END
