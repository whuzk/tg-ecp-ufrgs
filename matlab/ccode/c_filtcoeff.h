/*=========================================================================
 * c_filtcoeff.h
 * 
 *  Title: coefficients of Chebyshev and Butterworth filter for the
 *         elimination of mains interference (50 Hz or 60 Hz) and high
 *         frequency noise (>40 Hz).
 *  Author:     Diego Sogari
 *  Modified:   10/May/2014
 *
 *=======================================================================*/
#ifndef C_FILTCOEFF
#define C_FILTCOEFF

/*=========================================================================
 * Constants
 *=======================================================================*/
#define F50_COUNT   18
#define F60_COUNT   14
#define COEFF_COUNT 5

/*=========================================================================
 * Coefficients of 2nd order bandstop Chebyshev filter with Fc = 50 Hz
 * Sampling frequencies: i*50, i = 3,4,5,...,20
 *=======================================================================*/
static double chebb50[F50_COUNT][COEFF_COUNT] = {
    {0.959263793840377,   1.920211939798157,   2.879476473018586,   1.920211939798157,   0.959263793840377},
    {0.966521689593637,  -0.000000000000000,   1.933043379187275,  -0.000000000000000,   0.966521689593637},
    {0.970896599038088,  -1.200473317873331,   2.312877046848950,  -1.200473317873331,   0.970896599038088},
    {0.973821764366480,  -1.948070772845617,   2.921892584073222,  -1.948070772845617,   0.973821764366480},
    {0.975915409080945,  -2.434285461612747,   3.469827576251731,  -2.434285461612747,   0.975915409080945},
    {0.977487980574459,  -2.765094641234041,   3.910434372639535,  -2.765094641234041,   0.977487980574459},
    {0.978712485423220,  -2.999241397505042,   4.255201154987521,  -2.999241397505042,   0.978712485423220},
    {0.979692971331591,  -3.170603389822806,   4.524660510748602,  -3.170603389822806,   0.979692971331591},
    {0.980495772344351,  -3.299597437855944,   4.736970686820553,  -3.299597437855944,   0.980495772344351},
    {0.981165176744092,  -3.399042245119080,   4.906148697725432,  -3.399042245119080,   0.981165176744092},
    {0.981731883098379,  -3.477284103382074,   5.042589784156080,  -3.477284103382074,   0.981731883098379},
    {0.982217841405429,  -3.539933388490176,   5.153933945491846,  -3.539933388490176,   0.982217841405428},
    {0.982639162349158,  -3.590868182809498,   5.245814761768355,  -3.590868182809498,   0.982639162349158},
    {0.983007938005214,  -3.632835702085711,   5.322421945611917,  -3.632835702085711,   0.983007938005214},
    {0.983333421291929,  -3.667824637495369,   5.386904902260304,  -3.667824637495369,   0.983333421291929},
    {0.983622813050243,  -3.697302496893752,   5.441658168227678,  -3.697302496893752,   0.983622813050242},
    {0.983881801041169,  -3.722370899161813,   5.488523183229154,  -3.722370899161814,   0.983881801041169},
    {0.984114937534035,  -3.743869597281154,   5.528931736927929,  -3.743869597281154,   0.984114937534035}
};

static double cheba50[F50_COUNT][COEFF_COUNT] = {
    {1.000000000000000,   1.971586602749128,   2.911790521649109,   1.913307289499593,   0.941771807609906},
    {1.000000000000000,  -0.000000000000000,   1.954844233324163,  -0.000000000000000,   0.956009699833785},
    {1.000000000000000,  -1.225300525577104,   2.339284352166937,  -1.203447761918646,   0.964652652250254},
    {1.000000000000000,  -1.985403093996998,   2.955465652293106,  -1.955853644512659,   0.970457031198241},
    {1.000000000000000,  -2.478297386834696,   3.509814241506342,  -2.446648930707361,   0.974623998298189},
    {1.000000000000000,  -2.812840372094921,   3.955567805047302,  -2.781385484218079,   0.977760779416789},
    {1.000000000000000,  -3.049134479281939,   4.304357564007981,  -3.018807472881752,   0.980207381883923},
    {1.000000000000000,  -3.221743953091440,   4.576959154495125,  -3.192890540411604,   0.982169018242963},
    {1.000000000000000,  -3.351453621686167,   4.791744097134080,  -3.324156329241159,   0.983776881669439},
    {1.000000000000000,  -3.451289538967443,   4.962893531484692,  -3.425513059680259,   0.985118758737426},
    {1.000000000000000,  -3.529719777233525,   5.100924217067366,  -3.505378533978240,   0.986255609914229},
    {1.000000000000000,  -3.592428739697375,   5.213565374490657,  -3.569419030793135,   0.987231089327966},
    {1.000000000000000,  -3.643340806298313,   5.306516329773123,  -3.621556147072727,   0.988077279251870},
    {1.000000000000000,  -3.685232832131809,   5.384015667861730,  -3.664571081731824,   0.988818286142362},
    {1.000000000000000,  -3.720112845572537,   5.449249518361489,  -3.700479244664408,   0.989472574148772},
    {1.000000000000000,  -3.749461159596795,   5.504640332285893,  -3.730769324485759,   0.990054524644151},
    {1.000000000000000,  -3.774387967001210,   5.552051011503108,  -3.756559878452828,   0.990575505587552},
    {1.000000000000000,  -3.795738803229235,   5.592930030450957,  -3.778704324847151,   0.991044621510811}
};

static double chebd50[F50_COUNT] = {
    0.020221322489053,
    0.022781730722546,
    0.026399597464423,
    0.030419135518147,
    0.034637081731948,
    0.038968670178038,
    0.043371854587126,
    0.047823250171473,
    0.052308748853343,
    0.056819304747309,
    0.061348836389672,
    0.065893097123245,
    0.070449027787143,
    0.075014366589959,
    0.079587403943840,
    0.084166822877092,
    0.088751591971442,
    0.093340891649258
};

/*=========================================================================
 * Coefficients of 4th order lowpass Butterworth filter with Fc ~ 40 Hz
 * Sampling frequencies: i*50, i = 3,4,5,...,20
 *=======================================================================*/
static double buttb50[F50_COUNT][COEFF_COUNT] = {
    {0.115132028544437,   0.460528114177750,   0.690792171266625,   0.460528114177750,   0.115132028544437},
    {0.046582906636444,   0.186331626545775,   0.279497439818662,   0.186331626545775,   0.046582906636444},
    {0.022870207716291,   0.091480830865164,   0.137221246297746,   0.091480830865164,   0.022870207716291},
    {0.012634281782244,   0.050537127128977,   0.075805690693466,   0.050537127128977,   0.012634281782244},
    {0.007573673729384,   0.030294694917535,   0.045442042376303,   0.030294694917535,   0.007573673729384},
    {0.004824343357716,   0.019297373430865,   0.028946060146297,   0.019297373430865,   0.004824343357716},
    {0.003221793685119,   0.012887174740476,   0.019330762110713,   0.012887174740476,   0.003221793685119},
    {0.002234891698082,   0.008939566792329,   0.013409350188494,   0.008939566792329,   0.002234891698082},
    {0.001599557626992,   0.006398230507968,   0.009597345761952,   0.006398230507968,   0.001599557626992},
    {0.001175279549570,   0.004701118198282,   0.007051677297423,   0.004701118198282,   0.001175279549570},
    {0.000883064744212,   0.003532258976849,   0.005298388465274,   0.003532258976849,   0.000883064744212},
    {0.000676427630959,   0.002705710523837,   0.004058565785756,   0.002705710523837,   0.000676427630959},
    {0.000526933006904,   0.002107732027616,   0.003161598041423,   0.002107732027616,   0.000526933006904},
    {0.000416599204407,   0.001666396817626,   0.002499595226439,   0.001666396817626,   0.000416599204407},
    {0.000333721722732,   0.001334886890927,   0.002002330336390,   0.001334886890927,   0.000333721722732},
    {0.000270486149557,   0.001081944598227,   0.001622916897341,   0.001081944598227,   0.000270486149557},
    {0.000221556658605,   0.000886226634419,   0.001329339951629,   0.000886226634419,   0.000221556658605},
    {0.000183216023370,   0.000732864093478,   0.001099296140218,   0.000732864093478,   0.000183216023370}
};

static double butta50[F50_COUNT][COEFF_COUNT] = {
    {1.000000000000000,   0.260376974488161,   0.507465215378539,   0.055288349295379,   0.018981917548921},
    {1.000000000000000,  -0.782095198023338,   0.679978526916299,  -0.182675697753032,   0.030118875043169},
    {1.000000000000000,  -1.411983501196578,   1.122766080821220,  -0.408070951880240,   0.063211695716254},
    {1.000000000000000,  -1.835421688928288,   1.569754071766218,  -0.635723378097984,   0.103539503775964},
    {1.000000000000000,  -2.139938692798899,   1.968698669513034,  -0.853361038243604,   0.145779841199609},
    {1.000000000000000,  -2.369513007182038,   2.313988414415879,  -1.054665405878567,   0.187379492368185},
    {1.000000000000000,  -2.548772736075549,   2.611091129535982,  -1.237883637200263,   0.227113942701731},
    {1.000000000000000,  -2.692610987017437,   2.867399109111392,  -1.403484671368143,   0.264454816443505},
    {1.000000000000000,  -2.810570383161666,   3.089772909332597,  -1.552850946314607,   0.299241342175547},
    {1.000000000000000,  -2.909049493252733,   3.283996026668639,  -1.687645499468677,   0.331503438845899},
    {1.000000000000000,  -2.992499408722315,   3.454788121616806,  -1.809525101858776,   0.361365424871681},
    {1.000000000000000,  -3.064112180001714,   3.605962530284425,  -1.920021020979470,   0.388993512792107},
    {1.000000000000000,  -3.126236433358705,   3.740598954781649,  -2.020498979323391,   0.414567386010910},
    {1.000000000000000,  -3.180638548874721,   3.861194348994217,  -2.112155355110971,   0.438265142261981},
    {1.000000000000000,  -3.228672462827674,   3.969785202431404,  -2.196028886315913,   0.460255694275891},
    {1.000000000000000,  -3.271393338039506,   4.068043746665885,  -2.273017957826751,   0.480695327593280},
    {1.000000000000000,  -3.309635618293030,   4.157352813446632,  -2.343898842219977,   0.499726553604052},
    {1.000000000000000,  -3.344067837711877,   4.238863950884074,  -2.409342856586326,   0.517478199788043}
};

static double buttd50[F50_COUNT] = {
    1.185255627157917,
    1.808460866148262,
    2.388016376746198,
    2.947320984800637,
    3.495485083015091,
    4.036838554600198,
    4.573714898993687,
    5.107486854746922,
    5.639016200591772,
    6.168871876782752,
    6.697444917782047,
    7.225013089625836,
    7.751779206432059,
    8.277894854353058,
    8.803475628409958,
    9.328611229507803,
    9.853372338968903,
   10.377815411117203
};

/*=========================================================================
 * Coefficients of 2nd order bandstop Chebyshev filter with Fc = 60 Hz
 * Sampling frequencies: i*60, i = 3,4,5,...,16
 *=======================================================================*/
static double chebb60[F60_COUNT][COEFF_COUNT] = {
    {0.964097761390802,   1.929370843369529,   2.893468962965313,   1.929370843369529,   0.964097761390802},
    {0.970166382413556,  -0.000000000000000,   1.940332764827112,  -0.000000000000000,   0.970166382413556},
    {0.973821764366483,  -1.203973950108870,   2.319773554365469,  -1.203973950108870,   0.973821764366483},
    {0.976264695716705,  -1.952826816416280,   2.929091534786068,  -1.952826816416280,   0.976264695716705},
    {0.978012618664793,  -2.439396539363391,   3.477134311663941,  -2.439396539363391,   0.978012618664793},
    {0.979325197034834,  -2.770187279905455,   3.917636435934224,  -2.770187279905455,   0.979325197034834},
    {0.980347065671733,  -3.004161045601504,   4.262170798110865,  -3.004161045601504,   0.980347065671733},
    {0.981165176744081,  -3.175291312336781,   4.531335847460984,  -3.175291312336781,   0.981165176744081},
    {0.981834948951988,  -3.304038199172725,   4.743329664865065,  -3.304038199172725,   0.981834948951988},
    {0.982393374058915,  -3.403240058731440,   4.912191322290934,  -3.403240058731440,   0.982393374058915},
    {0.982866087968518,  -3.481251746726197,   5.048327524488828,  -3.481251746726196,   0.982866087968518},
    {0.983271417748266,  -3.543686878781159,   5.159383710400501,  -3.543686878781159,   0.983271417748266},
    {0.983622813050135,  -3.594424205164745,   5.250995554238549,  -3.594424205164745,   0.983622813050135},
    {0.983930367447167,  -3.636210393178489,   5.327353047791473,  -3.636210393178489,   0.983930367447167}
};

static double cheba60[F60_COUNT][COEFF_COUNT] = {
    {1.000000000000000,   1.976106412775101,   2.926256024114705,   1.927317396826315,   0.951240705001495},
    {1.000000000000000,  -0.000000000000000,   1.962394639137156,  -0.000000000000000,   0.963206879515664},
    {1.000000000000000,  -1.227046593459348,   2.346374412821837,  -1.208784029329181,   0.970457031198245},
    {1.000000000000000,  -1.987781263655173,   2.962827229904923,  -1.963097706773785,   0.975320202543883},
    {1.000000000000000,  -2.480857530335044,   3.517263987356528,  -2.454429309690372,   0.978808584121307},
    {1.000000000000000,  -2.815394548551956,   3.962898074436875,  -2.789134525056977,   0.981432950133077},
    {1.000000000000000,  -3.051604319557788,   4.311443376545143,  -3.026290862479055,   0.983478932725343},
    {1.000000000000000,  -3.224099290080563,   4.583740552964652,  -3.200019615619109,   0.985118758737430},
    {1.000000000000000,  -3.353686186347709,   4.798200507571794,  -3.330908130398327,   0.986462449855521},
    {1.000000000000000,  -3.453401065052123,   4.969026142413028,  -3.431894377601894,   0.987583581008425},
    {1.000000000000000,  -3.531716402561640,   5.106745465702225,  -3.511409081627701,   0.988533218201050},
    {1.000000000000000,  -3.594318304621304,   5.219092971687871,  -3.575123373189661,   0.989347914706920},
    {1.000000000000000,  -3.645131541090816,   5.311769949279581,  -3.626959810587686,   0.990054524644145},
    {1.000000000000000,  -3.686932730946118,   5.389015157141639,  -3.669698719264891,   0.990673219856634}
};

static double chebd60[F60_COUNT] = {
    0.016791742222555,
    0.018909535822730,
    0.021906477328454,
    0.025237750559759,
    0.028734245915009,
    0.032325397640336,
    0.035976192031983,
    0.039667151916624,
    0.043386525255853,
    0.047126776192506,
    0.050882838138257,
    0.054651172961156,
    0.058429231814358,
    0.062215130156192
};

/*=========================================================================
 * Coefficients of 4th order lowpass Butterworth filter with Fc ~ 40 Hz
 * Sampling frequencies: i*60, i = 3,4,5,...,16
 *=======================================================================*/
static double buttb60[F60_COUNT][COEFF_COUNT] = {
    {0.064920123876455,   0.259680495505819,   0.389520743258729,   0.259680495505819,   0.064920123876455},
    {0.026077721701092,   0.104310886804369,   0.156466330206554,   0.104310886804369,   0.026077721701092},
    {0.012634281782244,   0.050537127128977,   0.075805690693466,   0.050537127128977,   0.012634281782244},
    {0.006890401067214,   0.027561604268856,   0.041342406403284,   0.027561604268856,   0.006890401067214},
    {0.004084090581862,   0.016336362327449,   0.024504543491174,   0.016336362327449,   0.004084090581862},
    {0.002576434425323,   0.010305737701290,   0.015458606551936,   0.010305737701290,   0.002576434425323},
    {0.001706387732454,   0.006825550929817,   0.010238326394725,   0.006825550929817,   0.001706387732454},
    {0.001175279549570,   0.004701118198282,   0.007051677297423,   0.004701118198282,   0.001175279549570},
    {0.000835998871162,   0.003343995484647,   0.005015993226970,   0.003343995484647,   0.000835998871162},
    {0.000610957825789,   0.002443831303155,   0.003665746954733,   0.002443831303155,   0.000610957825789},
    {0.000456889466622,   0.001827557866489,   0.002741336799733,   0.001827557866489,   0.000456889466622},
    {0.000348518399541,   0.001394073598166,   0.002091110397249,   0.001394073598166,   0.000348518399541},
    {0.000270486149557,   0.001081944598227,   0.001622916897341,   0.001081944598227,   0.000270486149557},
    {0.000213138726975,   0.000852554907900,   0.001278832361850,   0.000852554907900,   0.000213138726975}
};

static double butta60[F60_COUNT][COEFF_COUNT] = {
    {1.000000000000000,  -0.434081628490660,   0.545644469464283,  -0.094204940981152,   0.021364082030805},
    {1.000000000000000,  -1.306605144101048,   1.030453835419574,  -0.362369044768858,   0.055763900667809},
    {1.000000000000000,  -1.835421688928288,   1.569754071766218,  -0.635723378097984,   0.103539503775964},
    {1.000000000000000,  -2.190866815260133,   2.041941424839012,  -0.895032246757243,   0.154204054253790},
    {1.000000000000000,  -2.446264193751164,   2.438200132788440,  -1.130128289436830,   0.203537799709350},
    {1.000000000000000,  -2.638627743891248,   2.769309786151488,  -1.339280761265205,   0.249821669810126},
    {1.000000000000000,  -2.788707999081582,   3.047719161665867,  -1.524197712136434,   0.292488753271416},
    {1.000000000000000,  -2.909049493252733,   3.283996026668639,  -1.687645499468677,   0.331503438845899},
    {1.000000000000000,  -3.007683635092584,   3.486484343171776,  -1.832488799851706,   0.367064073711100},
    {1.000000000000000,  -3.089990695845610,   3.661650221832193,  -1.961343652552220,   0.399459451778258},
    {1.000000000000000,  -3.159709994335829,   3.814501696447751,  -2.076481574331833,   0.429000103685866},
    {1.000000000000000,  -3.219520897144016,   3.948942541538915,  -2.179831812707754,   0.455986462705519},
    {1.000000000000000,  -3.271393338039506,   4.068043746665885,  -2.273017957826751,   0.480695327593280},
    {1.000000000000000,  -3.316807910624418,   4.174245550076572,  -2.357402780562258,   0.503375360741704}
};

static double buttd60[F60_COUNT] = {
    1.566738543382749,
    2.274162859032440,
    2.947320984800637,
    3.604206474885300,
    4.252034354248655,
    4.894287198027056,
    5.532859335169470,
    6.168871876782752,
    6.803031555102931,
    7.435806340485451,
    8.067518583112182,
    8.698397647770213,
    9.328611229507803,
    9.958284798181001
};

#endif