There are two types of files in this directory, both of which are the outputs from kinemetry- 

# Kinemetry_params
These files include the output parameters from the kinemetry run for each galaxy. The columns are:
KIN_PA, KIN_PA_e, alt_v_asym, alt_s_asym, K_asym, median(q)
<p>
- KIN_PA = the global kinematic position angle, this is the median PA while allowing the PA to vary
- KIN_PA_e = the error on the global position angle
- alt_v_asym = the sum of the higher order kinematic velocity moments, 
	sum = (k2 + k3 + k4 + k5)/4
        B_1_v = ABS(cf[i,2])
        elemental_v_asym[i]=sum/B_1_v
- alt_s_asym = the sum of the higher order kinematic moments, slightly different from v_asym, 
	sum = (k1 + k2 + k3 + k4 + k5)/5
        A_0_v = ABS(cf[i,0])
        elemental_s_asym[i]=sum/A_0_v
- K_asym = SQRT(alt_s_asym^2+alt_v_asym^2)
- median(q) = the median flattening of the best fit ellipses for each elliptical radius, cos i = q, this is a good way to find inclination

# Kinemetry_params_free
Same as above but allowing PA and q to vary at each elliptical radius; in this case, KIN_PA and KIN_PA_e are set to 0.

# vcirc_model and vcirc_model_free
These files contain the same format as the kinemetry_input_txt files but now have additional columns for vcirc and vkin.
<p>
xpos ypos (both of these are in pixels) vel_in vel_in_e vcirc vkin
<p>
- vcirc = circular velocity, the velocity field without any higher order terms
- vkin = adding in all higher order terms

Again, vcirc_model_free is for the fit where PA and q are allowed to vary with position.

