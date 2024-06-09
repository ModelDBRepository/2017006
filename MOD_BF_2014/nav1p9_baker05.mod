: nav1p9.mod is the NaV1.9 Na+ current from
: Baker 2005, parameter assignments and formula's from page 854

NEURON {
	SUFFIX nav1p9
	USEION na READ ena WRITE ina
	RANGE i, gbar, g
	RANGE tau_m,tau_h,minf,hinf
}

UNITS {
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.012  (S/cm2)
	: =100e-9/(100e-12*1e8) (S/cm2) : 100(nS)/100(um)^2

	A_am9 = 1.548 (/ms) : A for alpha m(9 etc ...)
	B_am9 = -11.01 (mV)
	C_am9 = -14.871 (mV)

	A_ah9 = 0.2574 (/ms) : A for alpha h
	B_ah9 = 63.264 (mV)
	C_ah9 = 3.7193 (mV)

	A_bm9 = 8.685 (/ms) : A for beta m
	B_bm9 = 112.4 (mV) : table has minus sign typo (Baker, personal comm.)
	C_bm9 = 22.9 (mV)

	A_bh9 = 0.53984 (/ms)   : A for beta h
	B_bh9 = 0.27853 (mV)
	C_bh9 = -9.0933 (mV)

	Q10m = 2.3 (1)	
	Q10h = 1.5 (1)
	Q10TEMP = 24 (degC)
}

ASSIGNED {
	v	(mV) : NEURON provides this
	ena	(mV)
	ina	(mA/cm2)
	tau_h	(ms)
	tau_m	(ms)
	minf
	hinf
	i	(mA/cm2)
	celsius (degC)
	qtm (1)
	qth (1)
	g (S/cm2)
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g  = gbar *m * h
	i = g * (v-ena)
	ina=i
}

INITIAL {
	: assume that equilibrium has been reached
	qtm = Q10m^((celsius-Q10TEMP)/10)
	qth = Q10h^((celsius-Q10TEMP)/10)
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
	h' = (hinf - h)/tau_h
}

FUNCTION alpham(Vm (mV)) (/ms) {
	alpham=A_am9/(1+exp((Vm+B_am9)/C_am9))
}

FUNCTION alphah(Vm (mV)) (/ms) {
	alphah=A_ah9/(1+exp((Vm+B_ah9)/C_ah9))
}

FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bm9/(1+exp((Vm+B_bm9)/C_bm9))
}

FUNCTION betah(Vm (mV)) (/ms) {
	betah=A_bh9/(1+exp((Vm+B_bh9)/C_bh9))
}

FUNCTION rates(Vm (mV)) (/ms) {
: Q10 only affects the tau, not inf state, Bin Feng 2014
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))/qtm
	minf = alpham(Vm) / (alpham(Vm) + betam(Vm))

	tau_h = 1.0 / (alphah(Vm) + betah(Vm))/qth
	hinf = alphah(Vm)/(alphah(Vm) + betah(Vm))
}
