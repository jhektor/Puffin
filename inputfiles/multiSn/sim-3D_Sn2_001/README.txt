	**	3/9	nuc_at_gvert diverged at time 24.77 dtmin reached...

	**	6/9 from advice on restart -> use checkpoints and reload all state vars.
			https://groups.google.com/forum/#!topic/moose-users/N3K6nTHJA5U

	** 	7/9 the following material blocks are altered using ? instead
			of GenericConstantMaterial
			diffusion_constants
			energy_chat
			energy_c_ab - not used
			energy_C
			energy_B
			energy_A
			ACMobility too???

	**	12/9	use 80x114 [x, y] elems for the mesh; square 16nm elems
						delta = 5x elem side_length

	**	13/9	wrong unit on ABC constants Vm was used to bring them into [Jm^-3] from [Jmol^-1]
						added postprocessing of the slip field TODO...
						grow time reduced from 30 to 8 cooldown 4 s -> low temp at 12 s

	**  14/9	added accum_L_2_slip_rate_gssi aux vars to store accumulated slip increments over time
						3D mesh - added layers in z-axis
