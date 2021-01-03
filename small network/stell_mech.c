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
 
#define nrn_init _nrn_init__stell_mech
#define _nrn_initial _nrn_initial__stell_mech
#define nrn_cur _nrn_cur__stell_mech
#define _nrn_current _nrn_current__stell_mech
#define nrn_jacob _nrn_jacob__stell_mech
#define nrn_state _nrn_state__stell_mech
#define _net_receive _net_receive__stell_mech 
#define _f_rates _f_rates__stell_mech 
#define rates rates__stell_mech 
#define states states__stell_mech 
 
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
#define gnap _p[0]
#define gnabar _p[1]
#define gkbar _p[2]
#define gks _p[3]
#define vhaks _p[4]
#define ghbar _p[5]
#define eh _p[6]
#define el _p[7]
#define gl _p[8]
#define gsynbar _p[9]
#define gna _p[10]
#define gk _p[11]
#define gh _p[12]
#define isyn _p[13]
#define il _p[14]
#define ih _p[15]
#define mna _p[16]
#define hna _p[17]
#define mnap _p[18]
#define nk _p[19]
#define mks _p[20]
#define mhf _p[21]
#define mhs _p[22]
#define msyn _p[23]
#define ena _p[24]
#define ek _p[25]
#define gsyn _p[26]
#define ina _p[27]
#define ik _p[28]
#define Dmna _p[29]
#define Dhna _p[30]
#define Dmnap _p[31]
#define Dnk _p[32]
#define Dmks _p[33]
#define Dmhf _p[34]
#define Dmhs _p[35]
#define Dmsyn _p[36]
#define _g _p[37]
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
 "setdata_stell_mech", _hoc_setdata,
 "nt_stell_mech", _hoc_nt,
 "rates_stell_mech", _hoc_rates,
 "vtrap_stell_mech", _hoc_vtrap,
 0, 0
};
#define nt nt_stell_mech
#define vtrap vtrap_stell_mech
 extern double nt( double );
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define asyn asyn_stell_mech
 double asyn = 1100;
#define ank ank_stell_mech
 double ank = 0;
#define amnap amnap_stell_mech
 double amnap = 0;
#define ahna ahna_stell_mech
 double ahna = 0;
#define amna amna_stell_mech
 double amna = 0;
#define bsyn bsyn_stell_mech
 double bsyn = 0.19;
#define bnk bnk_stell_mech
 double bnk = 0;
#define bmnap bmnap_stell_mech
 double bmnap = 0;
#define bhna bhna_stell_mech
 double bhna = 0;
#define bmna bmna_stell_mech
 double bmna = 0;
#define esyn esyn_stell_mech
 double esyn = 0;
#define flag flag_stell_mech
 double flag = 0;
#define mhstau mhstau_stell_mech
 double mhstau = 0;
#define mhsinf mhsinf_stell_mech
 double mhsinf = 0;
#define mhftau mhftau_stell_mech
 double mhftau = 0;
#define mhfinf mhfinf_stell_mech
 double mhfinf = 0;
#define mkstau mkstau_stell_mech
 double mkstau = 0;
#define mksinf mksinf_stell_mech
 double mksinf = 0;
#define usetable usetable_stell_mech
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_stell_mech", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "esyn_stell_mech", "mV",
 "asyn_stell_mech", "1/ms",
 "bsyn_stell_mech", "1/ms",
 "amna_stell_mech", "1/ms",
 "bmna_stell_mech", "1/ms",
 "ahna_stell_mech", "1/ms",
 "bhna_stell_mech", "1/ms",
 "amnap_stell_mech", "1/ms",
 "bmnap_stell_mech", "1/ms",
 "ank_stell_mech", "1/ms",
 "bnk_stell_mech", "1/ms",
 "mhstau_stell_mech", "ms",
 "mhftau_stell_mech", "ms",
 "mkstau_stell_mech", "ms",
 "gnap_stell_mech", "S/cm2",
 "gnabar_stell_mech", "S/cm2",
 "gkbar_stell_mech", "S/cm2",
 "gks_stell_mech", "S/cm2",
 "vhaks_stell_mech", "mV",
 "ghbar_stell_mech", "S/cm2",
 "eh_stell_mech", "mV",
 "el_stell_mech", "mV",
 "gl_stell_mech", "S/cm2",
 "gsynbar_stell_mech", "S/cm2",
 "gna_stell_mech", "S/cm2",
 "gk_stell_mech", "S/cm2",
 "gh_stell_mech", "S/cm2",
 "isyn_stell_mech", "mA/cm2",
 "il_stell_mech", "mA/cm2",
 "ih_stell_mech", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double hna0 = 0;
 static double msyn0 = 0;
 static double mhs0 = 0;
 static double mhf0 = 0;
 static double mks0 = 0;
 static double mnap0 = 0;
 static double mna0 = 0;
 static double nk0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "esyn_stell_mech", &esyn_stell_mech,
 "asyn_stell_mech", &asyn_stell_mech,
 "bsyn_stell_mech", &bsyn_stell_mech,
 "flag_stell_mech", &flag_stell_mech,
 "amna_stell_mech", &amna_stell_mech,
 "bmna_stell_mech", &bmna_stell_mech,
 "ahna_stell_mech", &ahna_stell_mech,
 "bhna_stell_mech", &bhna_stell_mech,
 "amnap_stell_mech", &amnap_stell_mech,
 "bmnap_stell_mech", &bmnap_stell_mech,
 "ank_stell_mech", &ank_stell_mech,
 "bnk_stell_mech", &bnk_stell_mech,
 "mksinf_stell_mech", &mksinf_stell_mech,
 "mhfinf_stell_mech", &mhfinf_stell_mech,
 "mhsinf_stell_mech", &mhsinf_stell_mech,
 "mhstau_stell_mech", &mhstau_stell_mech,
 "mhftau_stell_mech", &mhftau_stell_mech,
 "mkstau_stell_mech", &mkstau_stell_mech,
 "usetable_stell_mech", &usetable_stell_mech,
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
"stell_mech",
 "gnap_stell_mech",
 "gnabar_stell_mech",
 "gkbar_stell_mech",
 "gks_stell_mech",
 "vhaks_stell_mech",
 "ghbar_stell_mech",
 "eh_stell_mech",
 "el_stell_mech",
 "gl_stell_mech",
 "gsynbar_stell_mech",
 0,
 "gna_stell_mech",
 "gk_stell_mech",
 "gh_stell_mech",
 "isyn_stell_mech",
 "il_stell_mech",
 "ih_stell_mech",
 0,
 "mna_stell_mech",
 "hna_stell_mech",
 "mnap_stell_mech",
 "nk_stell_mech",
 "mks_stell_mech",
 "mhf_stell_mech",
 "mhs_stell_mech",
 "msyn_stell_mech",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 38, _prop);
 	/*initialize range parameters*/
 	gnap = 0.0005;
 	gnabar = 0.052;
 	gkbar = 0.011;
 	gks = 0;
 	vhaks = -35;
 	ghbar = 0.0015;
 	eh = -20;
 	el = -65;
 	gl = 0.0005;
 	gsynbar = 6e-006;
 	_prop->param = _p;
 	_prop->param_size = 38;
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

 void _stell_mech_reg() {
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
  hoc_register_prop_size(_mechtype, 38, 7);
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
 	ivoc_help("help ?1 stell_mech C:/Users/aniket/Desktop/NEURON 7.8 AMD64/stell/stell_mech.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_amna;
 static double *_t_bmna;
 static double *_t_ahna;
 static double *_t_bhna;
 static double *_t_amnap;
 static double *_t_bmnap;
 static double *_t_ank;
 static double *_t_bnk;
 static double *_t_mksinf;
 static double *_t_mkstau;
 static double *_t_mhfinf;
 static double *_t_mhftau;
 static double *_t_mhsinf;
 static double *_t_mhstau;
static int _reset;
static char *modelname = "Stell cells mechanism";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist1[8], _dlist1[8];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dmna = amna * ( 1.0 - mna ) - bmna * mna ;
   Dhna = ahna * ( 1.0 - hna ) - bhna * hna ;
   Dnk = ank * ( 1.0 - nk ) - bnk * nk ;
   Dmnap = amnap * ( 1.0 - mnap ) - bmnap * mnap ;
   Dmks = ( mksinf - mks ) / ( mkstau ) ;
   Dmhf = ( mhfinf - mhf ) / ( mhftau ) ;
   Dmhs = ( mhsinf - mhs ) / ( mhstau ) ;
   Dmsyn = asyn * nt ( _threadargscomma_ flag ) * ( 1.0 - msyn ) - bsyn * msyn ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dmna = Dmna  / (1. - dt*( ( amna )*( ( ( - 1.0 ) ) ) - ( bmna )*( 1.0 ) )) ;
 Dhna = Dhna  / (1. - dt*( ( ahna )*( ( ( - 1.0 ) ) ) - ( bhna )*( 1.0 ) )) ;
 Dnk = Dnk  / (1. - dt*( ( ank )*( ( ( - 1.0 ) ) ) - ( bnk )*( 1.0 ) )) ;
 Dmnap = Dmnap  / (1. - dt*( ( amnap )*( ( ( - 1.0 ) ) ) - ( bmnap )*( 1.0 ) )) ;
 Dmks = Dmks  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( mkstau ) )) ;
 Dmhf = Dmhf  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( mhftau ) )) ;
 Dmhs = Dmhs  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( mhstau ) )) ;
 Dmsyn = Dmsyn  / (1. - dt*( ( asyn * nt ( _threadargscomma_ flag ) )*( ( ( - 1.0 ) ) ) - ( bsyn )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    mna = mna + (1. - exp(dt*(( amna )*( ( ( - 1.0 ) ) ) - ( bmna )*( 1.0 ))))*(- ( ( amna )*( ( 1.0 ) ) ) / ( ( amna )*( ( ( - 1.0 ) ) ) - ( bmna )*( 1.0 ) ) - mna) ;
    hna = hna + (1. - exp(dt*(( ahna )*( ( ( - 1.0 ) ) ) - ( bhna )*( 1.0 ))))*(- ( ( ahna )*( ( 1.0 ) ) ) / ( ( ahna )*( ( ( - 1.0 ) ) ) - ( bhna )*( 1.0 ) ) - hna) ;
    nk = nk + (1. - exp(dt*(( ank )*( ( ( - 1.0 ) ) ) - ( bnk )*( 1.0 ))))*(- ( ( ank )*( ( 1.0 ) ) ) / ( ( ank )*( ( ( - 1.0 ) ) ) - ( bnk )*( 1.0 ) ) - nk) ;
    mnap = mnap + (1. - exp(dt*(( amnap )*( ( ( - 1.0 ) ) ) - ( bmnap )*( 1.0 ))))*(- ( ( amnap )*( ( 1.0 ) ) ) / ( ( amnap )*( ( ( - 1.0 ) ) ) - ( bmnap )*( 1.0 ) ) - mnap) ;
    mks = mks + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( mkstau ))))*(- ( ( ( mksinf ) ) / ( mkstau ) ) / ( ( ( ( - 1.0 ) ) ) / ( mkstau ) ) - mks) ;
    mhf = mhf + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( mhftau ))))*(- ( ( ( mhfinf ) ) / ( mhftau ) ) / ( ( ( ( - 1.0 ) ) ) / ( mhftau ) ) - mhf) ;
    mhs = mhs + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( mhstau ))))*(- ( ( ( mhsinf ) ) / ( mhstau ) ) / ( ( ( ( - 1.0 ) ) ) / ( mhstau ) ) - mhs) ;
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
    _t_amnap[_i] = amnap;
    _t_bmnap[_i] = bmnap;
    _t_ank[_i] = ank;
    _t_bnk[_i] = bnk;
    _t_mksinf[_i] = mksinf;
    _t_mkstau[_i] = mkstau;
    _t_mhfinf[_i] = mhfinf;
    _t_mhftau[_i] = mhftau;
    _t_mhsinf[_i] = mhsinf;
    _t_mhstau[_i] = mhstau;
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
  amnap = _xi;
  bmnap = _xi;
  ank = _xi;
  bnk = _xi;
  mksinf = _xi;
  mkstau = _xi;
  mhfinf = _xi;
  mhftau = _xi;
  mhsinf = _xi;
  mhstau = _xi;
  return;
 }
 if (_xi <= 0.) {
 amna = _t_amna[0];
 bmna = _t_bmna[0];
 ahna = _t_ahna[0];
 bhna = _t_bhna[0];
 amnap = _t_amnap[0];
 bmnap = _t_bmnap[0];
 ank = _t_ank[0];
 bnk = _t_bnk[0];
 mksinf = _t_mksinf[0];
 mkstau = _t_mkstau[0];
 mhfinf = _t_mhfinf[0];
 mhftau = _t_mhftau[0];
 mhsinf = _t_mhsinf[0];
 mhstau = _t_mhstau[0];
 return; }
 if (_xi >= 200.) {
 amna = _t_amna[200];
 bmna = _t_bmna[200];
 ahna = _t_ahna[200];
 bhna = _t_bhna[200];
 amnap = _t_amnap[200];
 bmnap = _t_bmnap[200];
 ank = _t_ank[200];
 bnk = _t_bnk[200];
 mksinf = _t_mksinf[200];
 mkstau = _t_mkstau[200];
 mhfinf = _t_mhfinf[200];
 mhftau = _t_mhftau[200];
 mhsinf = _t_mhsinf[200];
 mhstau = _t_mhstau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 amna = _t_amna[_i] + _theta*(_t_amna[_i+1] - _t_amna[_i]);
 bmna = _t_bmna[_i] + _theta*(_t_bmna[_i+1] - _t_bmna[_i]);
 ahna = _t_ahna[_i] + _theta*(_t_ahna[_i+1] - _t_ahna[_i]);
 bhna = _t_bhna[_i] + _theta*(_t_bhna[_i+1] - _t_bhna[_i]);
 amnap = _t_amnap[_i] + _theta*(_t_amnap[_i+1] - _t_amnap[_i]);
 bmnap = _t_bmnap[_i] + _theta*(_t_bmnap[_i+1] - _t_bmnap[_i]);
 ank = _t_ank[_i] + _theta*(_t_ank[_i+1] - _t_ank[_i]);
 bnk = _t_bnk[_i] + _theta*(_t_bnk[_i+1] - _t_bnk[_i]);
 mksinf = _t_mksinf[_i] + _theta*(_t_mksinf[_i+1] - _t_mksinf[_i]);
 mkstau = _t_mkstau[_i] + _theta*(_t_mkstau[_i+1] - _t_mkstau[_i]);
 mhfinf = _t_mhfinf[_i] + _theta*(_t_mhfinf[_i+1] - _t_mhfinf[_i]);
 mhftau = _t_mhftau[_i] + _theta*(_t_mhftau[_i+1] - _t_mhftau[_i]);
 mhsinf = _t_mhsinf[_i] + _theta*(_t_mhsinf[_i+1] - _t_mhsinf[_i]);
 mhstau = _t_mhstau[_i] + _theta*(_t_mhstau[_i+1] - _t_mhstau[_i]);
 }

 
static int  _f_rates (  double _lv ) {
    amna = ( .1 ) * vtrap ( _threadargscomma_ - ( _lv + 23.0 ) , 10.0 ) ;
   bmna = 4.0 * exp ( - ( _lv + 48.0 ) / 18.0 ) ;
   ahna = 0.07 * exp ( - ( _lv + 37.0 ) / 20.0 ) ;
   bhna = 1.0 / ( exp ( - 0.1 * ( _lv + 7.0 ) ) + 1.0 ) ;
   amnap = 1.0 / ( 0.15 * ( 1.0 + exp ( - ( _lv + 38.0 ) / 6.5 ) ) ) ;
   bmnap = ( exp ( - ( _lv + 38.0 ) / 6.5 ) ) / ( 0.15 * ( 1.0 + exp ( - ( _lv + 38.0 ) / 6.5 ) ) ) ;
   ank = 0.01 * vtrap ( _threadargscomma_ - ( _lv + 27.0 ) , 10.0 ) ;
   bnk = 0.125 * exp ( - ( _lv + 37.0 ) / 80.0 ) ;
   mksinf = 1.0 / ( 1.0 + exp ( - ( _lv - vhaks ) / 6.5 ) ) ;
   mkstau = 90.0 ;
   mhfinf = 1.0 / ( 1.0 + exp ( ( _lv + 79.2 ) / 9.78 ) ) ;
   mhftau = ( 0.51 / ( exp ( ( _lv - 1.7 ) / 10.0 ) ) ) + exp ( - ( _lv + 340.0 ) / 52.0 ) + 1.0 ;
   mhsinf = 1.0 / ( 1.0 + exp ( ( _lv + 71.3 ) / 7.9 ) ) ;
   mhstau = ( 5.6 / ( exp ( ( _lv - 1.7 ) / 14.0 ) ) ) + exp ( - ( _lv + 260.0 ) / 43.0 ) + 1.0 ;
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
 
static int _ode_count(int _type){ return 8;}
 
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
	for (_i=0; _i < 8; ++_i) {
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
  hna = hna0;
  msyn = msyn0;
  mhs = mhs0;
  mhf = mhf0;
  mks = mks0;
  mnap = mnap0;
  mna = mna0;
  nk = nk0;
 {
   rates ( _threadargscomma_ v ) ;
   mna = amna / ( amna + bmna ) ;
   hna = ahna / ( ahna + bhna ) ;
   nk = ank / ( ank + bnk ) ;
   mnap = amnap / ( amnap + bmnap ) ;
   msyn = asyn / ( asyn + bsyn ) ;
   mks = mksinf ;
   mhf = mhfinf ;
   mhs = mhsinf ;
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
   gna = ( gnabar * mna * mna * mna * hna ) + ( gnap * mnap ) ;
   ina = gna * ( v - ena ) ;
   gk = ( gkbar * nk * nk * nk * nk ) + ( gks * mks ) ;
   ik = gk * ( v - ek ) ;
   il = gl * ( v - el ) ;
   gsyn = gsynbar * msyn ;
   isyn = gsyn * ( v - esyn ) ;
   gh = ghbar * ( 0.65 * mhf + 0.35 * mhs ) ;
   ih = gh * ( v - eh ) ;
   }
 _current += ina;
 _current += ik;
 _current += il;
 _current += ih;
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
 if(error){fprintf(stderr,"at line 87 in file stell_mech.mod:\n\n"); nrn_complain(_p); abort_run(error);}
 }  }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(mna) - _p;  _dlist1[0] = &(Dmna) - _p;
 _slist1[1] = &(hna) - _p;  _dlist1[1] = &(Dhna) - _p;
 _slist1[2] = &(nk) - _p;  _dlist1[2] = &(Dnk) - _p;
 _slist1[3] = &(mnap) - _p;  _dlist1[3] = &(Dmnap) - _p;
 _slist1[4] = &(mks) - _p;  _dlist1[4] = &(Dmks) - _p;
 _slist1[5] = &(mhf) - _p;  _dlist1[5] = &(Dmhf) - _p;
 _slist1[6] = &(mhs) - _p;  _dlist1[6] = &(Dmhs) - _p;
 _slist1[7] = &(msyn) - _p;  _dlist1[7] = &(Dmsyn) - _p;
   _t_amna = makevector(201*sizeof(double));
   _t_bmna = makevector(201*sizeof(double));
   _t_ahna = makevector(201*sizeof(double));
   _t_bhna = makevector(201*sizeof(double));
   _t_amnap = makevector(201*sizeof(double));
   _t_bmnap = makevector(201*sizeof(double));
   _t_ank = makevector(201*sizeof(double));
   _t_bnk = makevector(201*sizeof(double));
   _t_mksinf = makevector(201*sizeof(double));
   _t_mkstau = makevector(201*sizeof(double));
   _t_mhfinf = makevector(201*sizeof(double));
   _t_mhftau = makevector(201*sizeof(double));
   _t_mhsinf = makevector(201*sizeof(double));
   _t_mhstau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "stell_mech.mod";
static const char* nmodl_file_text = 
  "TITLE Stell cells mechanism\n"
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
  "    SUFFIX stell_mech\n"
  "\n"
  "    USEION na READ ena WRITE ina \n"
  "    USEION k READ ek WRITE ik \n"
  "\n"
  "    NONSPECIFIC_CURRENT il \n"
  "    NONSPECIFIC_CURRENT ih \n"
  "\n"
  "    NONSPECIFIC_CURRENT isyn \n"
  "\n"
  "    RANGE gnabar,gkbar,gna,gk, gl,el, gksbar,gks,gnapbar,gnap,ghbar,gh,eh,vhaks,gsynbar\n"
  "    GLOBAL amna,bmna,ahna,bhna,amnap,bmnap,ank,bnk , mksinf, mkstau,mhfinf,mhftau,mhsinf,mhstau,flag\n"
  "    \n"
  "\n"
  "      \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    gnap = 0.0005 (S/cm2)\n"
  "    gnabar = 0.052 (S/cm2) \n"
  "\n"
  "    gkbar = 0.011 (S/cm2)\n"
  "    \n"
  "    gks = 0 (S/cm2)\n"
  "    vhaks = -35 (mV)\n"
  "\n"
  "    ghbar = 0.0015 (S/cm2)\n"
  "    eh = -20 (mV)\n"
  "\n"
  "    el = -65 (mV)\n"
  "    gl = 0.0005 (S/cm2)\n"
  "\n"
  "    gsynbar = 0.000006 (S/cm2)\n"
  "    esyn = 0 (mV)\n"
  "    asyn = 1100 (1/ms)\n"
  "    bsyn = 0.19 (1/ms)\n"
  "    flag = 0\n"
  "\n"
  "    \n"
  "    \n"
  "\n"
  "\n"
  "\n"
  "\n"
  "}   \n"
  "\n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "        ena (mV)\n"
  "        ek (mV)\n"
  "\n"
  "	    gna (S/cm2)\n"
  "	    gk (S/cm2)\n"
  "        gh (S/cm2)\n"
  "        gsyn (S/cm2)\n"
  "\n"
  "        isyn (mA/cm2)\n"
  "        ina (mA/cm2)\n"
  "        ik (mA/cm2)\n"
  "        il (mA/cm2)\n"
  "        ih (mA/cm2)\n"
  "\n"
  "    amna (1/ms) bmna (1/ms) ahna (1/ms) bhna (1/ms) amnap (1/ms) bmnap (1/ms) ank (1/ms) bnk (1/ms) \n"
  "    mksinf  mhfinf mhsinf      \n"
  "    mhstau (ms)  mhftau (ms)   mkstau (ms)\n"
  "} \n"
  "\n"
  "\n"
  "STATE {\n"
  "    mna hna mnap nk mks mhf mhs msyn\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT { \n"
  "    SOLVE states METHOD cnexp\n"
  "\n"
  "    gna = (gnabar*mna*mna*mna*hna) + (gnap*mnap)\n"
  "	ina = gna*(v - ena)\n"
  "\n"
  "    gk = (gkbar*nk*nk*nk*nk) + (gks * mks)\n"
  "	ik = gk*(v - ek)  \n"
  "\n"
  "    il = gl*(v - el)\n"
  "\n"
  "    gsyn = gsynbar*msyn\n"
  "    isyn = gsyn*(v-esyn)\n"
  "\n"
  "    gh = ghbar*(0.65*mhf+0.35*mhs)\n"
  "    ih = gh*(v-eh)\n"
  "    \n"
  "\n"
  "\n"
  "\n"
  "}\n"
  "INITIAL {\n"
  "    rates(v)\n"
  "    mna = amna/(amna+bmna)\n"
  "    hna = ahna/(ahna+bhna)\n"
  "    nk = ank/(ank+bnk)\n"
  "    mnap = amnap/(amnap+bmnap)\n"
  "    msyn = asyn/(asyn+bsyn)\n"
  "    mks = mksinf\n"
  "    mhf = mhfinf\n"
  "    mhs = mhsinf\n"
  "    \n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates(v)\n"
  "\n"
  "    mna' = amna*(1-mna) - bmna*mna\n"
  "    hna' = ahna*(1-hna) - bhna*hna\n"
  "    nk' = ank*(1-nk) - bnk*nk\n"
  "    mnap' = amnap*(1-mnap) - bmnap*mnap\n"
  "    mks' = (mksinf-mks)/(mkstau)\n"
  "    mhf' = (mhfinf-mhf)/(mhftau)\n"
  "    mhs' = (mhsinf-mhs)/(mhstau)\n"
  "    msyn' = asyn*nt(flag)*(1-msyn) - bsyn*msyn    \n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v (mV)) {\n"
  "\n"
  "    TABLE amna,bmna,ahna,bhna,amnap,bmnap,ank,bnk,mksinf, mkstau,mhfinf,mhftau,mhsinf,mhstau FROM -100 TO 100 WITH 200 \n"
  "    UNITSOFF\n"
  "\n"
  "    amna = (.1)*vtrap(-(v+23),10)\n"
  "    bmna = 4*exp(-(v+48)/18)\n"
  "\n"
  "\n"
  "    ahna = 0.07*exp(-(v+37)/20)\n"
  "    bhna = 1/(exp(-0.1*(v+7))+1)\n"
  "\n"
  "\n"
  "    amnap = 1/(0.15*(1+exp(-(v+38)/6.5)))\n"
  "    bmnap = (exp(-(v+38)/6.5))/(0.15*(1+exp(-(v+38)/6.5)))\n"
  "\n"
  " \n"
  "    ank = 0.01*vtrap(-(v+27),10)\n"
  "    bnk = 0.125 * exp(-(v+37)/80)\n"
  "\n"
  "\n"
  "    mksinf = 1/(1+exp(-(v-vhaks)/6.5))\n"
  "    mkstau = 90\n"
  "\n"
  "    mhfinf = 1/(1+exp((v+79.2)/9.78))\n"
  "    mhftau = (0.51 / (exp((v-1.7)/10))) + exp(-(v+340)/52) + 1\n"
  "\n"
  "    mhsinf = 1/(1+exp((v+71.3)/7.9))\n"
  "    mhstau = (5.6/ (exp((v-1.7)/14))) + exp(-(v+260)/43) + 1\n"
  "    \n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator..\n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{\n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
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
