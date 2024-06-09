: slow-adapting ms channel
: gating probability from Haselwandter and Phillips 2013
: PLOS Computational Biology 9(5): e1003055 (2013).
: relative Ina and Ik were calculated based upon GHK equation assuming equal permeability
NEURON {
	SUFFIX ms_v1
	USEION na READ nai, nao, ena WRITE ina
	USEION k READ ki, ko, ek WRITE ik
	RANGE p_m, p_m_inf, tau_m, tension0, slope_m, A_m
	RANGE m1, m2, tau_t
	RANGE gbar

	RANGE lamda, dlamda, tension
	RANGE Pna2k, im

}

UNITS {

	(S) = (siemens) 
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.01 (S/cm2)
	m1=0.2 (1)
	m2=0.6 (ms)
	tau_t = 30 (ms)
	slope_m =2.09 (1)
	A_m = 300 (ms)

: ::::::::::::::::::: the variables below all need to be played during simulation 
	lamda = 1 (1)
	dlamda = 0 (1/ms)
: ::::::::::::::::::::::::::::::::::

	tension0 = 10.45 (1)

	Q10=3 (1)	
	Q10TEMP = 24 (degC)
}

ASSIGNED {
	v       (mV)
	i	  (mA/cm2)	

	ena	  (mV)
	nai	  (mM)
	nao	  (mM)
	ina	  (mA/cm2)
	ek	  (mV)
	ki	  (mM)
	ko	  (mM)
	ik	  (mA/cm2)

	im	  (mA/cm2)

	Pna2k    (1)
	p_m_inf  (1)
	tau_m	  (ms)

}

STATE { tension (1) 
	p_m (1) }

BREAKPOINT {
	SOLVE state METHOD cnexp

::::::::: portions of Ina and Ik is determined by their relative driving force
	Pna2k = (v-ena)/(v-ek)
	if (Pna2k < 0 ) {Pna2k = - Pna2k}

	ik =gbar*p_m*(1/(1+Pna2k))*(v - ek)
	ina =gbar*p_m*(Pna2k/(1+Pna2k))*(v - ena)
	im=ik+ina
}

DERIVATIVE state {
	p_m_inf = 1/(1+exp((tension0-tension)/slope_m))
	tau_m = A_m/(exp((tension -tension0)/slope_m/2)+exp((tension0 -tension)/slope_m/2))

	p_m' = (p_m_inf - p_m)/tau_m
	tension' = (m1*lamda - m1 +  m2*dlamda - tension)/tau_t
}

INITIAL { 
	tension = 0
	p_m = 1/(1+exp(tension0/slope_m))
}
