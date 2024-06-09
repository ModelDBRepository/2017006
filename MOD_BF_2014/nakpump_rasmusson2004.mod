: sodium-potassium ATPase adapted from Randall Rasmusson 2004

NEURON {
	SUFFIX nakpump
	USEION na READ nai, nao WRITE ina
	USEION k READ ko WRITE ik
	RANGE Imax, Km_nai, Km_ko
	RANGE inak, fnak
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolts)
}

PARAMETER {
	Imax = .0667 (mA/cm2) <0, 1e6>: BR model trsducer zone 0.2 mA/cm2
	Km_nai = 21 (mM) <0, 1e6>
	Km_ko = 1.5 (mM) <0, 1e6>

	Q10=1.16 (1)	
	Q10TEMP = 24 (degC)

	RdivF = 0.0861728 (mV/degC)
}

ASSIGNED {
	ko (mM)
	nai (mM)
	nao (mM)
	v       (mV)

	celsius (degC)
	qt (1)
	fnak (1)
	inak (mA/cm2)
	ik (mA/cm2)
	ina (mA/cm2)
}

BREAKPOINT {
inak = finak(v)
ina = 3.0*inak
ik = -2.0*inak

}

INITIAL { 
	qt = Q10^((celsius-Q10TEMP)/10)
}

FUNCTION finak(Vm (mV)) (mA/cm2) {
LOCAL delta1, x1
UNITSOFF
	x1 = Vm/(RdivF*(273.15+celsius))
	delta1 = (exp(nao/67.3)-1.0)/7.0
	fnak = 1.0/(1.0+0.1245*exp(-0.1*x1)+0.0365*delta1*exp(-1.0*x1))
	finak = Imax * fnak *(1.0/(1+(Km_nai/nai)^4))*ko/(ko+Km_ko)*qt
UNITSON
}


