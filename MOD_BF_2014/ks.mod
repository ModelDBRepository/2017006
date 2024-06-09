: ks.mod is the sustained K+ current from
: value optimized from the GUI channel builder 2014

NEURON {
	SUFFIX ks
	USEION k READ ek WRITE ik
	RANGE gbar : parameters that can be modified by hoc
	RANGE i : calculated variables accessible by hoc
}

UNITS {
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.36

	A_an = 0.1(/ms)  
	K_an = 0.08 (/mV)
	D_an = -52 (mV)

	A_bn = 0.125 (/ms)
	K_bn = -0.0125 (/mV)
	D_bn = -60 (mV)


	tau_noff = 0.0 (ms)

	Q10=1.4 (1)	
	Q10TEMP = 24 (degC)

}

ASSIGNED {
	v	(mV) : NEURON provides this
	ek	(mV)
	ik	(mA/cm2)
	i	(mA/cm2)
	tau_n	(ms)
	ninf

	celsius (degC)
	qt (1)
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gbar*n^4* (v - ek)
	ik=i
}

INITIAL {
	: assume that equilibrium has been reached
	qt = Q10^((celsius-Q10TEMP)/10)
	rates(v)
	n = ninf
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/tau_n
}

FUNCTION alphan(Vm (mV)) (/ms) {
	UNITSOFF
	if (Vm - D_an != 0) {
		alphan=A_an*K_an*(Vm - D_an)/(1-exp(-K_an*(Vm - D_an)))
	} else {
		alphan=A_an
	}	
	UNITSON
}

FUNCTION betan(Vm (mV)) (/ms) {
	UNITSOFF
		betan=A_bn*exp(K_bn*(Vm - D_bn))
	UNITSON

}

FUNCTION rates(Vm (mV)) (/ms) {
	UNITSOFF
qt = Q10^((celsius-Q10TEMP)/10)
	tau_n = (1.0 / (alphan(Vm) + betan(Vm)) + tau_noff)/qt
	ninf = alphan(Vm) / (alphan(Vm) + betan(Vm))
	UNITSON
}
