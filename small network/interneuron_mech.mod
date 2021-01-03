TITLE interneuron mech



UNITS {
    (mV)=(millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX interneuron_mech


    USEION na READ ena WRITE ina 
    USEION k READ ek WRITE ik 

    NONSPECIFIC_CURRENT il 
     

    NONSPECIFIC_CURRENT isyn 

    RANGE gnabar,gkbar,gna,gk, gl,el,gsynbar
    GLOBAL amna,bmna,ahna,bhna,ank,bnk,minf,hinf,ninf,mtau,htau,ntau,flag
    

      
}

PARAMETER {
        gnabar = 0.035 (S/cm2)	
        gkbar = 0.009 (S/cm2)

        gl = 0.0003 (S/cm2)	
        el = -54.3 (mV)

		gsynbar = 0.000006 (S/cm2)
    	esyn = 0 (mV)
    	asyn = 3.33 (1/ms)
    	bsyn = 0.11 (1/ms)
    	flag = 0

}

STATE {
        mna hna nk msyn m h n q10 sum
}
 
ASSIGNED {
        v (mV)
        ena (mV)
        celsius (degC)
        ek (mV)

	    gna (S/cm2)
	    gk (S/cm2)
	    gsyn (S/cm2)

	    isyn (mA/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)

    minf hinf ninf
	mtau (ms) htau (ms) ntau (ms)

	amna (1/ms) bmna (1/ms) ahna (1/ms) bhna (1/ms) ank (1/ms) bnk (1/ms)

}


BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	    ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)      
        il = gl*(v - el)
        gsyn = gsynbar*msyn
    	isyn = gsyn*(v-esyn)





}

INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	msyn = asyn/(asyn+bsyn)
}


DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
        msyn' = asyn*nt(flag)*(1-msyn) - bsyn*msyn

        
}


PROCEDURE rates(v (mV)) {
     
    TABLE amna,bmna,ahna,bhna,ank,bnk,minf, mtau, hinf, htau, ninf, ntau FROM -100 TO 100 WITH 200
    
    UNITSOFF
    q10 = 3^((celsius - 6.3)/10)
            :"m" sodium activation system
    amna = (.1)*vtrap(-(v+23),10)
    bmna = 4*exp(-(v+60)/18)
    sum = amna + bmna
    mtau = 1/(q10*sum)
    minf = amna/sum


    ahna = 0.07*exp(-(v+58)/20)
    bhna = 1/(exp(-0.1*(v+28))+1)
    sum = ahna + bhna
    htau = 1/(q10*sum)
    hinf = ahna/sum


 
    ank = 0.01*vtrap(-(v+27),10)
    bnk = 0.125 * exp(-(v+44)/80)
    sum = ank + bnk
    ntau = 1/(q10*sum)
    ninf = ank/sum
    


}






FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}



FUNCTION nt(flag) {
    if (flag == 0){
        nt=0
    }
    else {
        nt= 0.001
    }
}
UNITSON
