TITLE Stell cells mechanism



UNITS {
    (mV)=(millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX stell_mech

    USEION na READ ena WRITE ina 
    USEION k READ ek WRITE ik 

    NONSPECIFIC_CURRENT il 
    NONSPECIFIC_CURRENT ih 

    NONSPECIFIC_CURRENT isyn 

    RANGE gnabar,gkbar,gna,gk, gl,el, gksbar,gks,gnapbar,gnap,ghbar,gh,eh,vhaks,gsynbar
    GLOBAL amna,bmna,ahna,bhna,amnap,bmnap,ank,bnk , mksinf, mkstau,mhfinf,mhftau,mhsinf,mhstau,flag
    

      
}

PARAMETER {
    gnap = 0.0005 (S/cm2)
    gnabar = 0.052 (S/cm2) 

    gkbar = 0.011 (S/cm2)
    
    gks = 0 (S/cm2)
    vhaks = -35 (mV)

    ghbar = 0.0015 (S/cm2)
    eh = -20 (mV)

    el = -65 (mV)
    gl = 0.0005 (S/cm2)

    gsynbar = 0.000006 (S/cm2)
    esyn = 0 (mV)
    asyn = 1100 (1/ms)
    bsyn = 0.19 (1/ms)
    flag = 0

    
    




}   

ASSIGNED {
        v (mV)
        ena (mV)
        ek (mV)

	    gna (S/cm2)
	    gk (S/cm2)
        gh (S/cm2)
        gsyn (S/cm2)

        isyn (mA/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        ih (mA/cm2)

    amna (1/ms) bmna (1/ms) ahna (1/ms) bhna (1/ms) amnap (1/ms) bmnap (1/ms) ank (1/ms) bnk (1/ms) 
    mksinf  mhfinf mhsinf      
    mhstau (ms)  mhftau (ms)   mkstau (ms)
} 


STATE {
    mna hna mnap nk mks mhf mhs msyn
}


BREAKPOINT { 
    SOLVE states METHOD cnexp

    gna = (gnabar*mna*mna*mna*hna) + (gnap*mnap)
	ina = gna*(v - ena)

    gk = (gkbar*nk*nk*nk*nk) + (gks * mks)
	ik = gk*(v - ek)  

    il = gl*(v - el)

    gsyn = gsynbar*msyn
    isyn = gsyn*(v-esyn)

    gh = ghbar*(0.65*mhf+0.35*mhs)
    ih = gh*(v-eh)
    



}
INITIAL {
    rates(v)
    mna = amna/(amna+bmna)
    hna = ahna/(ahna+bhna)
    nk = ank/(ank+bnk)
    mnap = amnap/(amnap+bmnap)
    msyn = asyn/(asyn+bsyn)
    mks = mksinf
    mhf = mhfinf
    mhs = mhsinf
    

}

DERIVATIVE states {
    rates(v)

    mna' = amna*(1-mna) - bmna*mna
    hna' = ahna*(1-hna) - bhna*hna
    nk' = ank*(1-nk) - bnk*nk
    mnap' = amnap*(1-mnap) - bmnap*mnap
    mks' = (mksinf-mks)/(mkstau)
    mhf' = (mhfinf-mhf)/(mhftau)
    mhs' = (mhsinf-mhs)/(mhstau)
    msyn' = asyn*nt(flag)*(1-msyn) - bsyn*msyn    


}

PROCEDURE rates(v (mV)) {

    TABLE amna,bmna,ahna,bhna,amnap,bmnap,ank,bnk,mksinf, mkstau,mhfinf,mhftau,mhsinf,mhstau FROM -100 TO 100 WITH 200 
    UNITSOFF

    amna = (.1)*vtrap(-(v+23),10)
    bmna = 4*exp(-(v+48)/18)


    ahna = 0.07*exp(-(v+37)/20)
    bhna = 1/(exp(-0.1*(v+7))+1)


    amnap = 1/(0.15*(1+exp(-(v+38)/6.5)))
    bmnap = (exp(-(v+38)/6.5))/(0.15*(1+exp(-(v+38)/6.5)))

 
    ank = 0.01*vtrap(-(v+27),10)
    bnk = 0.125 * exp(-(v+37)/80)


    mksinf = 1/(1+exp(-(v-vhaks)/6.5))
    mkstau = 90

    mhfinf = 1/(1+exp((v+79.2)/9.78))
    mhftau = (0.51 / (exp((v-1.7)/10))) + exp(-(v+340)/52) + 1

    mhsinf = 1/(1+exp((v+71.3)/7.9))
    mhstau = (5.6/ (exp((v-1.7)/14))) + exp(-(v+260)/43) + 1
    


}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator..
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
