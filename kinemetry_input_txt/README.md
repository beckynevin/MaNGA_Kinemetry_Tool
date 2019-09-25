The columns of the kinemetry input file are:
<p>
#	xbin	ybin	velbin	velbin_e	sigbin	sigbin_e
<p>
- # = idk why but kinemetry requires that the first column is the number of the entry.
<p>
- x/ybin = x/y position in pixels relative to the kinematic center
<p>
- velbin, velbin_e = stellar velocity and error at that position
<p>
- sigbin, sigbin_e = stellar velocity dispersion and error at that position
<p>
* Note - Voronoi or similar binning schemes are very popular to deal with low S/N IFS data. Kinemetry can still function with binned data, simply feed it the same value for velbin, velbin_e, sigbin, and sigbin_e for each spaxel in a given bin.
