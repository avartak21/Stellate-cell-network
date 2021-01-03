/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__interneuron_mech
#define _nrn_initial _nrn_initial__interneuron_mech
#define nrn_cur _nrn_cur__interneuron_mech
#define _nrn_current _nrn_current__interneuron_mech
#define nrn_jacob _nrn_jacob__interneuron_mech
#define nrn_state _nrn_state__interneuron_mech
#define _net_receive _net_receive__interneuron_mech 
#define _f_rates _f_rates__interneuron_mech 
#define rates rates__interneuron_mech 
#define states states__interneuron_mech 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnabar _p[0]
#define gkbar _p[1]
#define gl _p[2]
#define el _p[3]
#define gsynbar _p[4]
#define gna _p[5]
#define gk _p[6]
#define isyn _p[7]
#define il _p[8]
#define mna _p[9]
#define hna _p[10]
#define nk _p[11]
#define msyn _p[12]
#define m _p[13]
#define h _p[14]
#define n _p[15]
#define q10 _p[16]
#define sum _p[17]
#define Dmna _p[18]
#define Dhna _p[19]
#define Dnk _p[20]
#define Dmsyn _p[21]
#define Dm _p[22]
#define Dh _p[23]
#define Dn _p[24]
#define Dq10 _p[25]
#define Dsum _p[26]
#define ena _p[27]
#define ek _p[28]
#define gsyn _p[29]
#define ina _p[30]
#define ik _p[31]
#define _g _p[32]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_nt(void);
 static void _hoc_rates(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_interneuron_mech", _hoc_setdata,
 "nt_interneuron_mech", _hoc_nt,
 "rates_interneuron_mech", _hoc_rates,
 "vtrap_interneuron_mech", _hoc_vtrap,
 0, 0
};
#define nt nt_interneuron_mech
#define vtrap vtrap_interneuron_mech
 extern double nt( double );
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define asyn asyn_interneuron_mech
 double asyn = 3.33;
#define ank ank_interneuron_mech
 double ank = 0;
#define ahna ahna_interneuron_mech
 double ahna = 0;
#define amna amna_interneuron_mech
 double amna = 0;
#define bsyn bsyn_interneuron_mech
 double bsyn = 0.11;
#define bnk bnk_interneuron_mech
 double bnk = 0;
#define bhna bhna_interneuron_mech
 double bhna = 0;
#define bmna bmna_interneuron_mech
 double bmna = 0;
#define esyn esyn_interneuron_mech
 double esyn = 0;
#define flag flag_interneuron_mech
 double flag = 0;
#define htau htau_interneuron_mech
 double htau = 0;
#define hinf hinf_interneuron_mech
 double hinf = 0;
#define mtau mtau_interneuron_mech
 double mtau = 0;
#define minf minf_interneuron_mech
 double minf = 0;
#define ntau ntau_interneuron_mech
 double ntau = 0;
#define ninf ninf_interneuron_mech
 double ninf = 0;
#define usetable usetable_interneuron_mech
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_interneuron_mech", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "esyn_interneuron_mech", "mV",
 "asyn_interneuron_mech", "1/ms",
 "bsyn_interneuron_mech", "1/ms",
 "mtau_interneuron_mech", "ms",
 "htau_interneuron_mech", "ms",
 "ntau_interneuron_mech", "ms",
 "amna_interneuron_mech", "1/ms",
 "bmna_interneuron_mech", "1/ms",
 "ahna_interneuron_mech", "1/ms",
 "bhna_interneuron_mech", "1/ms",
 "ank_interneuron_mech", "1/ms",
 "bnk_interneuron_mech", "1/ms",
 "gnabar_interneuron_mech", "S/cm2",
 "gkbar_interneuron_mech", "S/cm2",
 "gl_interneuron_mech", "S/cm2",
 "el_interneuron_mech", "mV",
 "gsynbar_interneuron_mech", "S/cm2",
 "gna_interneuron_mech", "S/cm2",
 "gk_interneuron_mech", "S/cm2",
 "isyn_interneuron_mech", "mA/cm2",
 "il_interneuron_mech", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double hna0 = 0;
 static double m0 = 0;
 static double msyn0 = 0;
 static double mna0 = 0;
 static double n0 = 0;
 static double nk0 = 0;
 static double q100 = 0;
 static double sum0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "esyn_interneuron_mech", &esyn_interneuron_mech,
 "asyn_interneuron_mech", &asyn_interneuron_mech,
 "bsyn_interneuron_mech", &bsyn_interneuron_mech,
 "flag_interneuron_mech", &flag_interneuron_mech,
 "minf_interneuron_mech", &minf_interneuron_mech,
 "hinf_interneuron_mech", &hinf_interneuron_mech,
 "ninf_interneuron_mech", &ninf_interneuron_mech,
 "mtau_interneuron_mech", &mtau_interneuron_mech,
 "htau_interneuron_mech", &htau_interneuron_mech,
 "ntau_interneuron_mech", &ntau_interneuron_mech,
 "amna_interneuron_mech", &amna_interneuron_mech,
 "bmna_interneuron_mech", &bmna_interneuron_mech,
 "ahna_interneuron_mech", &ahna_interneuron_mech,
 "bhna_interneuron_mech", &bhna_interneuron_mech,
 "ank_interneuron_mech", &ank_interneuron_mech,
 "bnk_interneuron_mech", &bnk_interneuron_mech,
 "usetable_interneuron_mech", &usetable_interneuron_mech,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"interneuron_mech",
 "gnabar_interneuron_mech",
 "gkbar_interneuron_mech",
 "gl_interneuron_mech",
 "el_interneuron_mech",
 "gsynbar_interneuron_mech",
 0,
 "gna_interneuron_mech",
 "gk_interneuron_mech",
 "isyn_interneuron_mech",
 "il_interneuron_mech",
 0,
 "mna_interneuron_mech",
 "hna_interneuron_mech",
 "nk_interneuron_mech",
 "msyn_interneuron_mech",
 "m_interneuron_mech",
 "h_interneuron_mech",
 "n_interneuron_mech",
 "q10_interneuron_mech",
 "sum_interneuron_mech",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 33, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.035;
 	gkbar = 0.009;
 	gl = 0.0003;
 	el = -54.3;
 	gsynbar = 6e-006;
 	_prop->param = _p;
 	_prop->param_size = 33;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _interneuron_mech_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 33, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 interneuron_mech C:/Users/aniket/Desktop/NEURON 7.8 AMD64/interneuron/interneuron_mech.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_amna;
 static double *_t_bmna;
 static double *_t_ahna;
 static double *_t_bhna;
 static double *_t_ank;
 static double *_t_bnk;
 static double *_t_minf;
 static double *_t_mtau;
 static double *_t_hinf;
 static double *_t_htau;
 static double *_t_ninf;
 static double *_t_ntau;
static int _reset;
static char *modelname = "interneuron mech";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist1[4], _dlist1[4];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dn = ( ninf - n ) / ntau ;
   Dmsyn = asyn * nt ( _threadargscomma_ flag ) * ( 1.0 - msyn ) - bsyn * msyn ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
 Dmsyn = Dmsyn  / (1. - dt*( ( asyn * nt ( _threadargscomma_ flag ) )*( ( ( - 1.0 ) ) ) - ( bsyn )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ntau)))*(- ( ( ( ninf ) ) / ntau ) / ( ( ( ( - 1.0 ) ) ) / ntau ) - n) ;
    msyn = msyn + (1. - exp(dt*(( asyn * nt ( _threadargscomma_ flag ) )*( ( ( - 1.0 ) ) ) - ( bsyn )*( 1.0 ))))*(- ( ( ( asyn )*( nt ( _threadargscomma_ flag ) ) )*( ( 1.0 ) ) ) / ( ( ( asyn )*( nt ( _threadargscomma_ flag ) ) )*( ( ( - 1.0 ) ) ) - ( bsyn )*( 1.0 ) ) - msyn) ;
   }
  return 0;
}
 static double _mfac_rates, _tmin_rates;
 static void _check_rates();
 static void _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_x);
    _t_amna[_i] = amna;
    _t_bmna[_i] = bmna;
    _t_ahna[_i] = ahna;
    _t_bhna[_i] = bhna;
    _t_ank[_i] = ank;
    _t_bnk[_i] = bnk;
    _t_minf[_i] = minf;
    _t_mtau[_i] = mtau;
    _t_hinf[_i] = hinf;
    _t_htau[_i] = htau;
    _t_ninf[_i] = ninf;
    _t_ntau[_i] = ntau;
   }
  }
 }

 static int rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return 0;
 }

 static void _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  amna = _xi;
  bmna = _xi;
  ahna = _xi;
  bhna = _xi;
  ank = _xi;
  bnk = _xi;
  minf = _xi;
  mtau = _xi;
  hinf = _xi;
  htau = _xi;
  ninf = _xi;
  ntau = _xi;
  return;
 }
 if (_xi <= 0.) {
 amna = _t_amna[0];
 bmna = _t_bmna[0];
 ahna = _t_ahna[0];
 bhna = _t_bhna[0];
 ank = _t_ank[0];
 bnk = _t_bnk[0];
 minf = _t_minf[0];
 mtau = _t_mtau[0];
 hinf = _t_hinf[0];
 htau = _t_htau[0];
 ninf = _t_ninf[0];
 ntau = _t_ntau[0];
 return; }
 if (_xi >= 200.) {
 amna = _t_amna[200];
 bmna = _t_bmna[200];
 ahna = _t_ahna[200];
 bhna = _t_bhna[200];
 ank = _t_ank[200];
 bnk = _t_bnk[200];
 minf = _t_minf[200];
 mtau = _t_mtau[200];
 hinf = _t_hinf[200];
 htau = _t_htau[200];
 ninf = _t_ninf[200];
 ntau = _t_ntau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 amna = _t_amna[_i] + _theta*(_t_amna[_i+1] - _t_amna[_i]);
 bmna = _t_bmna[_i] + _theta*(_t_bmna[_i+1] - _t_bmna[_i]);
 ahna = _t_ahna[_i] + _theta*(_t_ahna[_i+1] - _t_ahna[_i]);
 bhna = _t_bhna[_i] + _theta*(_t_bhna[_i+1] - _t_bhna[_i]);
 ank = _t_ank[_i] + _theta*(_t_ank[_i+1] - _t_ank[_i]);
 bnk = _t_bnk[_i] + _theta*(_t_bnk[_i+1] - _t_bnk[_i]);
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 }

 
static int  _f_rates (  double _lv ) {
    q10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   amna = ( .1 ) * vtrap ( _threadargscomma_ - ( _lv + 23.0 ) , 10.0 ) ;
   bmna = 4.0 * exp ( - ( _lv + 60.0 ) / 18.0 ) ;
   sum = amna + bmna ;
   mtau = 1.0 / ( q10 * sum ) ;
   minf = amna / sum ;
   ahna = 0.07 * exp ( - ( _lv + 58.0 ) / 20.0 ) ;
   bhna = 1.0 / ( exp ( - 0.1 * ( _lv + 28.0 ) ) + 1.0 ) ;
   sum = ahna + bhna ;
   htau = 1.0 / ( q10 * sum ) ;
   hinf = ahna / sum ;
   ank = 0.01 * vtrap ( _threadargscomma_ - ( _lv + 27.0 ) , 10.0 ) ;
   bnk = 0.125 * exp ( - ( _lv + 44.0 ) / 80.0 ) ;
   sum = ank + bnk ;
   ntau = 1.0 / ( q10 * sum ) ;
   ninf = ank / sum ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
    _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap (  double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   _r =  vtrap (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double nt (  double _lflag ) {
   double _lnt;
 if ( _lflag  == 0.0 ) {
     _lnt = 0.0 ;
     }
   else {
     _lnt = 0.001 ;
     }
   
return _lnt;
 }
 
static void _hoc_nt(void) {
  double _r;
   _r =  nt (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
     _ode_spec1 ();
   }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  hna = hna0;
  m = m0;
  msyn = msyn0;
  mna = mna0;
  n = n0;
  nk = nk0;
  q10 = q100;
  sum = sum0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   n = ninf ;
   msyn = asyn / ( asyn + bsyn ) ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
  ek = _ion_ek;
 initmodel();
  }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gna = gnabar * m * m * m * h ;
   ina = gna * ( v - ena ) ;
   gk = gkbar * n * n * n * n ;
   ik = gk * ( v - ek ) ;
   il = gl * ( v - el ) ;
   gsyn = gsynbar * msyn ;
   isyn = gsyn * ( v - esyn ) ;
   }
 _current += ina;
 _current += ik;
 _current += il;
 _current += isyn;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 73 in file interneuron_mech.mod:\n        SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }  }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(n) - _p;  _dlist1[2] = &(Dn) - _p;
 _slist1[3] = &(msyn) - _p;  _dlist1[3] = &(Dmsyn) - _p;
   _t_amna = makevector(201*sizeof(double));
   _t_bmna = makevector(201*sizeof(double));
   _t_ahna = makevector(201*sizeof(double));
   _t_bhna = makevector(201*sizeof(double));
   _t_ank = makevector(201*sizeof(double));
   _t_bnk = makevector(201*sizeof(double));
   _t_minf = makevector(201*sizeof(double));
   _t_mtau = makevector(201*sizeof(double));
   _t_hinf = makevector(201*sizeof(double));
   _t_htau = makevector(201*sizeof(double));
   _t_ninf = makevector(201*sizeof(double));
   _t_ntau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "interneuron_mech.mod";
static const char* nmodl_file_text = 
  "TITLE interneuron mech\n"
  "\n"
  "\n"
  "\n"
  "UNITS {\n"
  "    (mV)=(millivolt)\n"
  "    (S) = (siemens)\n"
  "    (mA) = (milliamp)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX interneuron_mech\n"
  "\n"
  "\n"
  "    USEION na READ ena WRITE ina \n"
  "    USEION k READ ek WRITE ik \n"
  "\n"
  "    NONSPECIFIC_CURRENT il \n"
  "     \n"
  "\n"
  "    NONSPECIFIC_CURRENT isyn \n"
  "\n"
  "    RANGE gnabar,gkbar,gna,gk, gl,el,gsynbar\n"
  "    GLOBAL amna,bmna,ahna,bhna,ank,bnk,minf,hinf,ninf,mtau,htau,ntau,flag\n"
  "    \n"
  "\n"
  "      \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "        gnabar = 0.035 (S/cm2)	\n"
  "        gkbar = 0.009 (S/cm2)\n"
  "\n"
  "        gl = 0.0003 (S/cm2)	\n"
  "        el = -54.3 (mV)\n"
  "\n"
  "		gsynbar = 0.000006 (S/cm2)\n"
  "    	esyn = 0 (mV)\n"
  "    	asyn = 3.33 (1/ms)\n"
  "    	bsyn = 0.11 (1/ms)\n"
  "    	flag = 0\n"
  "\n"
  "}\n"
  "\n"
  "STATE {\n"
  "        mna hna nk msyn m h n q10 sum\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "        ena (mV)\n"
  "        celsius (degC)\n"
  "        ek (mV)\n"
  "\n"
  "	    gna (S/cm2)\n"
  "	    gk (S/cm2)\n"
  "	    gsyn (S/cm2)\n"
  "\n"
  "	    isyn (mA/cm2)\n"
  "        ina (mA/cm2)\n"
  "        ik (mA/cm2)\n"
  "        il (mA/cm2)\n"
  "\n"
  "    minf hinf ninf\n"
  "	mtau (ms) htau (ms) ntau (ms)\n"
  "\n"
  "	amna (1/ms) bmna (1/ms) ahna (1/ms) bhna (1/ms) ank (1/ms) bnk (1/ms)\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "        gna = gnabar*m*m*m*h\n"
  "	    ina = gna*(v - ena)\n"
  "        gk = gkbar*n*n*n*n\n"
  "	    ik = gk*(v - ek)      \n"
  "        il = gl*(v - el)\n"
  "        gsyn = gsynbar*msyn\n"
  "    	isyn = gsyn*(v-esyn)\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m = minf\n"
  "	h = hinf\n"
  "	n = ninf\n"
  "	msyn = asyn/(asyn+bsyn)\n"
  "}\n"
  "\n"
  "\n"
  "DERIVATIVE states {  \n"
  "        rates(v)\n"
  "        m' =  (minf-m)/mtau\n"
  "        h' = (hinf-h)/htau\n"
  "        n' = (ninf-n)/ntau\n"
  "        msyn' = asyn*nt(flag)*(1-msyn) - bsyn*msyn\n"
  "\n"
  "        \n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE rates(v (mV)) {\n"
  "     \n"
  "    TABLE amna,bmna,ahna,bhna,ank,bnk,minf, mtau, hinf, htau, ninf, ntau FROM -100 TO 100 WITH 200\n"
  "    \n"
  "    UNITSOFF\n"
  "    q10 = 3^((celsius - 6.3)/10)\n"
  "            :\"m\" sodium activation system\n"
  "    amna = (.1)*vtrap(-(v+23),10)\n"
  "    bmna = 4*exp(-(v+60)/18)\n"
  "    sum = amna + bmna\n"
  "    mtau = 1/(q10*sum)\n"
  "    minf = amna/sum\n"
  "\n"
  "\n"
  "    ahna = 0.07*exp(-(v+58)/20)\n"
  "    bhna = 1/(exp(-0.1*(v+28))+1)\n"
  "    sum = ahna + bhna\n"
  "    htau = 1/(q10*sum)\n"
  "    hinf = ahna/sum\n"
  "\n"
  "\n"
  " \n"
  "    ank = 0.01*vtrap(-(v+27),10)\n"
  "    bnk = 0.125 * exp(-(v+44)/80)\n"
  "    sum = ank + bnk\n"
  "    ntau = 1/(q10*sum)\n"
  "    ninf = ank/sum\n"
  "    \n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.\n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{\n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "FUNCTION nt(flag) {\n"
  "    if (flag == 0){\n"
  "        nt=0\n"
  "    }\n"
  "    else {\n"
  "        nt= 0.001\n"
  "    }\n"
  "}\n"
  "UNITSON\n"
  ;
#endif
