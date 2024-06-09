: kf.mod is the slow inactivating K+ current
: from schild 1994 & Bin Feng 2008

NEURON {
	SUFFIX kf_d
	USEION k READ ek WRITE ik
	RANGE gbar
	RANGE i
	RANGE tau_x, xinf, tau_y, yinf 
	RANGE B_xinf, B_yinf
}

UNITS {
	(S) = (siemens)	
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.003 : =30e-9/(100e-12*1e8) (S/cm2) : 30(nS)/100(um)^2

	A_xinf = 1.0   : for xinf
	B_xinf = -14.59 (mV)   : original -14.59 (mV)
	C_xinf = -15 (mV)

	A_xtau = 5.0 (ms) : for tau_x
	B_xtau = -65 (mV)
	C_xtau = 0.022 (/mV)
	D_xtau = 3.5 (ms)

	A_yinf = 1.0   : for yinf
	B_yinf = -48.0 (mV)
	C_yinf = 7.0 (mV)

	D_ytau = 1800 (ms) : for tau_y

	Q10=1.93 (1)	
	Q10TEMP = 24 (degC)

}

ASSIGNED {
	v	(mV) : NEURON provides this
	ek	(mV)
	ik	(mA/cm2)
	i	(mA/cm2)
	tau_x	(ms)
	tau_y	(ms)
	xinf
	yinf
	celsius (degC)
	qt (1)
}

STATE {x1 y1}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gbar*x1*x1*x1*y1*(v-ek)
	ik = i
}

INITIAL {
	: assume that equilibrium has been reached
	qt = Q10^((celsius-Q10TEMP)/10)
	rates(v)
	x1 = xinf
	y1 = yinf
}

DERIVATIVE states {
	rates(v)
	x1' = (xinf - x1)/tau_x
	y1' = (yinf - y1)/tau_y
}



FUNCTION rates(Vm (mV)) (/ms) {
UNITSOFF

qt = Q10^((celsius-Q10TEMP)/10)
	xinf = A_xinf/(1+exp((Vm - B_xinf)/C_xinf))
	yinf = A_yinf/(1+exp((Vm - B_yinf)/C_yinf))
	tau_x = (A_xtau * exp (-(Vm - B_xtau)^2*(C_xtau)^2) + D_xtau)/qt
	tau_y = D_ytau/qt

UNITSON
}
