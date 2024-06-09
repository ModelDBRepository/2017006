: passive leak current

NEURON {
	SUFFIX nakpas
	USEION k READ ek WRITE ik
	USEION na READ ena WRITE ina

	RANGE gbar, epas, i
	RANGE gna_p,gk_p : passive relative peameability Na: K is 0.1:1
}

UNITS {
	(S) = (siemens)	
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 1e-3 (S/cm2)
	epas = -60 (mV)
}

ASSIGNED {
	v	(mV) : NEURON provides this
	gna_p	(S/cm2)
	gk_p	(S/cm2)
	ek	(mV)
	ena	(mV)
	i	(mA/cm2)
	ik	(mA/cm2)
	ina	(mA/cm2)

}

BREAKPOINT {
	gna_p=gbar*(epas-ek)/(ena-ek)
	gk_p=gbar*(ena-epas)/(ena-ek)
	ina = gna_p*(v-ena)
	ik = gk_p*(v-ek)
	i = ina + ik
}

