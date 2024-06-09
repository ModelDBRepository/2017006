:Tracking intracellular Na and K concentration
: transmembrane Na+ and K+ flux
:Longitudinal diffusion of Na+ and K+ to adjacent compartments (no buffering)
:(equivalent modified euler with standard method and
:equivalent to diagonalized linear solver with CVODE )
: Adapted from Fleidervish 2010

NEURON {
	SUFFIX ioni
	USEION na READ ina WRITE nai
	USEION k READ ik WRITE ki
	RANGE D
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = .6 (um2/ms)
}

ASSIGNED {
	ina (milliamp/cm2)
	ik (milliamp/cm2)
	diam (um)
}

STATE {
	nai (mM)
	ki  (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	COMPARTMENT PI*diam*diam/4 {nai ki}
	LONGITUDINAL_DIFFUSION D*PI*diam*diam/4 {nai ki}
	~ nai << (-ina/(FARADAY)*PI*diam*(1e4))
	~ ki << (-ik/(FARADAY)*PI*diam*(1e4))
}
