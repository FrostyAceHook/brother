#include "nozzle.h"

#include "assertion.h"
#include "maths.h"


// Nozzle diverging and exit angles extrapolated from data in nasa paper:
// https://ntrs.nasa.gov/citations/19770009165

#define NZL_IN_AEAT(AEAT) (within((AEAT), NZL_min_AEAT, NZL_max_AEAT))
#define NZL_IN_NLF(NLF) (within((NLF), NZL_min_NLF, NZL_max_NLF))

#define NZL_ASSERT_IN_AEATNLF(AEAT, NLF) do {                           \
        assert(NZL_IN_AEAT((AEAT)), "AEAT=%f, NLF=%f", (AEAT), (NLF));  \
        assert(NZL_IN_NLF((NLF)), "AEAT=%f, NLF=%f", (AEAT), (NLF));    \
    } while (0)

f64 nzl_theta_div(f64 AEAT, f64 NLF) {
    NZL_ASSERT_IN_AEATNLF(AEAT, NLF);
    f64 x = log2(AEAT);
    f64 y = NLF;
    f64 A;
    f64 B;
    f64 C;
    f64 D;
    f64 E;
    f64 F;
    f64 G;
    if (NLF >= 0.9) {
        A = -73.576668144;
        B = 2.0140420776;
        C = 89.1204667856;
        D = -278.203452801;
        E = 24.1209238324;
        F = 344.177751142;
        G = -24.3867292652;
    } else if (NLF >= 0.8) {
        A = -36.7821906186;
        B = 5.44659097911;
        C = 56.9174532511;
        D = -159.494910997;
        E = 19.8640728454;
        F = 253.581035759;
        G = -16.8130123402;
    } else if (NLF >= 0.7) {
        A = 31.6379252304;
        B = 55.8899435791;
        C = 59.7342379608;
        D = -0.628213679815;
        E = 31.902201812;
        F = 492.781259534;
        G = 32.8799713668;
    } else {
        A = -11.8198005202;
        B = 45.2314009746;
        C = 102.053219128;
        D = -161.375493686;
        E = 38.7116965704;
        F = 629.775421116;
        G = 8.3754483057;
    }
    return (A + B*x + C*y + x*y)
         / (D + E*x + F*y + G*x*y);
}

f64 nzl_theta_exit(f64 AEAT, f64 NLF) {
    NZL_ASSERT_IN_AEATNLF(AEAT, NLF);
    f64 x = log2(AEAT);
    f64 y = NLF;
    f64 A;
    f64 B;
    f64 C;
    f64 D;
    f64 E;
    f64 F;
    f64 G;
    if (NLF >= 0.9) {
        A = -77.1832953717;
        B = 4.89002019523;
        C = 111.684097893;
        D = 108.741300367;
        E = -457.430200479;
        F = -128.559316968;
        G = 611.611528301;
    } else if (NLF >= 0.8) {
        A = 85.2862055798;
        B = -0.0581506237535;
        C = -90.9930605548;
        D = 189.894422987;
        E = 59.9863440541;
        F = -212.118540681;
        G = -51.6239165065;
    } else if (NLF >= 0.7) {
        A = -136.251317086;
        B = 2.56249489274;
        C = 241.089075404;
        D = -239.854862081;
        E = -270.48519269;
        F = 414.264012684;
        G = 443.983190791;
    } else {
        A = -68.1116220938;
        B = 1.75535969394;
        C = 132.256411039;
        D = -160.591697636;
        E = -60.2164539955;
        F = 283.313697674;
        G = 129.355138767;
    }
    return (A + B*x + C*y + x*y)
         / (D + E*x + F*y + G*x*y);
}
