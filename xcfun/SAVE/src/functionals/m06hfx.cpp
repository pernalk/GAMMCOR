#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M06-HF exchange functional. to be used with HF exchange factor of 1.00 

template<class num>
static num energy (const densvars<num> &d)
{
   using pw91_like_x_internal::chi2;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::fw;
   using m0xy_metagga_xc_internal::h;
   using m0xy_metagga_xc_internal::alpha_x;

     const parameter param_a[12] =
       {  1.179732e-01, -1.066708e+00, -1.462405e-01,  7.481848e+00,  3.776679e+00,
         -4.436118e+01, -1.830962e+01,  1.003903e+02,  3.864360e+01, -9.806018e+01,
         -2.557716e+01,  3.590404e+01 };
     const parameter param_d[6] =
       { -1.179732e-01, -2.500000e-03, -1.180056e-02,  0.000000e+00,  0.000000e+00,
          0.000000e+00 };

   num chia2 = chi2(d.a, d.gaa);
   num chib2 = chi2(d.b, d.gbb);

   return (  (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua) 
              + lsda_x(d.a)*h(param_d, alpha_x, chia2, zet(d.a, d.taua)))
            +
             (  pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) 
              + lsda_x(d.b)*h(param_d, alpha_x, chib2, zet(d.b, d.taub)))
          );
}

void setup_m06hfx(functional &f)
{
  f.describe(XC_M06HFX, XC_MGGA,
	     "M06-HF exchange",
             "M06-HF Meta-Hybrid Exchange Functional\n"
             "Y Zhao and D. G. Truhlar, Theor. Chem. Account 120, 215 (2008)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);
  const double d[] =
    {1., .8, 1., 1., 1.,  0.165,   0.1050};
  const double ref[] =
    {  0.19601260,  0.87809662,  0.90304189,  0.00284062, 
       0.00000000, 0.00398338,  -2.7579263,  -3.42255772 };
  f.add_test(XC_VARS_AB,1,d,ref,1e-5);
}
