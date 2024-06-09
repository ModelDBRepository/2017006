: kf_a.mod is the fast inactivating K+ current
: from schild 1994 & Bin Feng 2008, A-type K current

NEURON {
	SUFFIX kf_a
	USEION k READ ek WRITE ik
	RANGE gbar
	RANGE i
	RANGE tau_p, pinf, tau_q, qinf 
}

UNITS {
	(S) = (siemens)	
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.006 : =30e-9/(100e-12*1e8) (S/cm2) : 30(nS)/100(um)^2

	A_pinf = 1.0   : for pinf
	B_pinf = -15.79 (mV) 		: original 5mV
	C_pinf = -10 (mV)

	A_ptau = 5.0 (ms) : for tau_p
	B_ptau = -65 (mV)
	C_ptau = 0.022 (/mV)
	D_ptau = 1.5 (ms)

	A_qinf = 1.0   : for qinf
	B_qinf = -58.0 (mV)
	C_qinf = 7.0 (mV)

	A_qtau = 45.0 (ms) : for tau_q
	B_qtau = -10 (mV)
	C_qtau = 0.0035 (/mV)
	D_qtau = 10.5 (ms)

	Q10=1.93 (1)	
	Q10TEMP = 24 (degC)

}

ASSIGNED {
	v	(mV) : NEURON provides this
	ek	(mV)
	ik	(mA/cm2)
	i	(mA/cm2)
	tau_p	(ms)
	tau_q	(ms)
	pinf
	qinf
	celsius (degC)
	qt (1)
}

STATE {p q}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gbar*p*p*p*q*(v-ek)
	ik = i
}

INITIAL {
	: assume that equilibrium has been reached
	qt = Q10^((celsius-Q10TEMP)/10)
	rates(v)
	p = pinf
	q = qinf
}

DERIVATIVE states {
	rates(v)
	p' = (pinf - p)/tau_p
	q' = (qinf - q)/tau_q
}



FUNCTION rates(Vm (mV)) (/ms) {
UNITSOFF
qt = Q10^((celsius-Q10TEMP)/10)
	pinf = A_pinf/(1+exp((Vm - B_pinf)/C_pinf))
	qinf = A_qinf/(1+exp((Vm - B_qinf)/C_qinf))
	tau_p = (A_ptau * exp (-(Vm - B_ptau)^2*(C_ptau)^2) + D_ptau)/qt
	tau_q = (A_qtau * exp (-(Vm - B_qtau)^2*(C_qtau)^2) + D_qtau)/qt
UNITSON
}
