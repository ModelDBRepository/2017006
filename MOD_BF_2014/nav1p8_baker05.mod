: nav1p8.mod is the NaV1.8 Na+ current from
: Baker 2005, parameter assignments and formula's from page 854

NEURON {
	SUFFIX nav1p8
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
	gbar = 0.012 (S/cm2)
	: =220e-9/(100e-12*1e8) (S/cm2) : 220(nS)/100(um)^2

	A_am8 = 3.83 (/ms) : A for alpha m(8 etc ...)
	B_am8 = 2.58 (mV)
	C_am8 = -11.47 (mV)

	A_ah8 = 0.013536 (/ms) : A for alpha h
	B_ah8 = 105 (mV)
	C_ah8 = 46.33 (mV)

	A_bm8 = 6.894 (/ms) : A for beta m
	B_bm8 = 61.2 (mV)
	C_bm8 = 19.8 (mV)

	A_bh8 = 0.61714 (/ms)   : A for beta h
	B_bh8 = -21.8 (mV)
	C_bh8 = -11.998 (mV)

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
	g  (S/cm2)
	
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3 * h
	i = g * (v-ena)
	ina = i
}

INITIAL {
	qtm = Q10m^((celsius-Q10TEMP)/10)
	qth = Q10h^((celsius-Q10TEMP)/10)
	rates(v)
	m=minf
	h=hinf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
	h' = (hinf - h)/tau_h
}

FUNCTION alpham(Vm (mV)) (/ms) {
	alpham=A_am8/(1+exp((Vm+B_am8)/C_am8))
}

FUNCTION alphah(Vm (mV)) (/ms) {
	alphah=A_ah8*exp(-(Vm+B_ah8)/C_ah8)
}

FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bm8/(1+exp((Vm+B_bm8)/C_bm8))
}

FUNCTION betah(Vm (mV)) (/ms) {
	betah=A_bh8/(1+exp((Vm+B_bh8)/C_bh8))
}

FUNCTION rates(Vm (mV)) (/ms) {
: Q10 only affects the tau, not inf state, Bin Feng 2014
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))/qtm
	minf = alpham(Vm) / (alpham(Vm) + betam(Vm))

	tau_h = 1.0 / (alphah(Vm) + betah(Vm))/qth
	hinf = alphah(Vm) / (alphah(Vm) + betah(Vm))
}
