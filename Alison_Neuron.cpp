#include <Alison_Neuron.h>
#include <math.h>

// Model Functions

static inline double
m_inf(double V)
{
    //1.0/(1.0+exp(-(VV-thetam)/sigmam))
    return 1.0/(1.0+exp(-(V-(-38.0))/5.0));
}

static inline double
h_inf(double V)
{
    //1.0/(1.0+exp(-(VV-thetah)/sigmah))
    return 1.0/(1.0+exp(-(V-(-50.0))/-3.0));
}

static inline double
n_inf(double V)
{
    //1.0/(1.0+exp(-(VV-thetan)/sigman))
    return 1.0/(1.0+exp(-(V - (-40.0))/5.0));
}

static inline double
tau_m(double V)
{
    //0.05+0.5*(GAMMAF(V,thetatma,sigmatma)*GAMMAF(V,thetatmb,sigmatmb))
    return 0.05+0.5*((1.0/(1.0+exp(-(V - (-20.0))/-10.0)))*(1.0/(1.0+exp(-(V - (-60.0))/3.0))));
}

static inline double
tau_h(double V)
{
    //1+8*GAMMAF(V,thetath,sigmath)
    return 1.0+8.0*(1.0/(1.0+exp(-(V - (-45.0))/-3.0)));
}

static inline double
tau_n(double V)
{
    //1.2+8*(GAMMAF(V,thetatna,sigmatna)*GAMMAF(V,thetatnb,sigmatnb))
    return 1.2+8.0*((1.0/(1.0+exp(-(V - (-50.0))/-10.0)))*(1.0/(1.0+exp(-(V - (-70.0))/10.0))));
}

static inline double
tau_ns(double V)
{
    //1.2+200*(GAMMAF(V,thetatna,sigmatna)*GAMMAF(V,thetatnb,sigmatnb))
    return 1.2+200*((1.0/(1.0+exp(-(V - (-50.0))/-10.0)))*(1.0/(1.0+exp(-(V - (-70.0))/10.0))));
}

static inline double
a_inf(double V)
{
    //ainf(v)=1/(1+exp(-(v-thetaa)/sigmaa))
    return 1.0/(1.0+exp(-(V - (-50.0))/20.0));
}

static inline double
b_inf(double V)
{
    //abinf(v)=1/(1+exp(-(v-thetab)/sigmab))
    return 1.0/(1.0+exp(-(V - (-70.0))/-6.0));
}

static inline double
tau_b(double V,double taub1)
{
    //taub(v)=taub0 + (taub1-taub0)/(1+exp(-(v-tb)/sb))
    return 10 + (taub1-10)/(1.0+exp(-(V-(-80))/10));
}

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
    return new Alison_Neuron();
}

static DefaultGUIModel::variable_t vars[] =
{
    { "Vm", "Membrane Potential", DefaultGUIModel::OUTPUT, },
    { "GA_INPUT", "Input from GA_Calc", DefaultGUIModel::INPUT, },
    { "IA-ACTIVATION", "Input from ia-activate", DefaultGUIModel::INPUT, },
    { "Iapp (nA)", "Applied Current",
        DefaultGUIModel::PARAMETER, },
    { "V0 (mV)", "Initial membrane potential (mV)",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "cm (nF)", "Specific membrane capacitance",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_Na_max (uS)", "Maximum Na+ conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "E_Na (mV)", "Sodium reversal potential",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_K_max (uS)",
        "Maximum delayed rectifier conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_Ks_max (uS)",
        "Maximum delayed rectifier conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },        
    { "E_K (mV)", "K+ reversal potential", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::DOUBLE, },
    { "G_A_max (uS)",
        "Maximum transient A-type K+ conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "G_L (uS)", "Maximum leak conductance density",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "E_L (mV)", "Leak reversal potential",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "Rate (Hz)", "Rate of integration (Hz)", DefaultGUIModel::PARAMETER
        | DefaultGUIModel::UINTEGER, },
    { "taua", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "taub1", "",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
    { "m", "", DefaultGUIModel::STATE, },
    { "h", "", DefaultGUIModel::STATE, },
    { "n", "", DefaultGUIModel::STATE, },
    { "ns", "", DefaultGUIModel::STATE, },
    { "a", "", DefaultGUIModel::STATE, },
    { "b", "", DefaultGUIModel::STATE, },
    { "IA", "", DefaultGUIModel::STATE, },
    { "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, };

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

/*
* Macros for making the code below a little bit cleaner.
*/

#define V (y[0])
#define m (y[1])
#define h (y[2])
#define n (y[3])
#define ns (y[4])
#define b (y[5])
#define a (y[6])
#define dV (dydt[0])
#define dm (dydt[1])
#define dh (dydt[2])
#define dn (dydt[3])
#define dns (dydt[4])
#define db (dydt[5])
#define da (dydt[6])
#define G_Na (G_Na_max*m*m*m*h)
#define G_K  (G_K_max*n*n*n*n)
#define G_Ks (G_Ks_max*ns*ns*ns*ns)
#define G_A  (G_A_max*a*a*a*b)

Alison_Neuron::Alison_Neuron(void) : DefaultGUIModel("Alison_Neuron", ::vars, ::num_vars) {
    setWhatsThis(
    "<p><b>Alison_Neuron:</b><br>This module simulates a neuron with an A Current.</p>");
    createGUI(vars, num_vars);
    initParameters();
    update( INIT );
    refresh();
    resizeMe();
}

Alison_Neuron::~Alison_Neuron(void) {}

void Alison_Neuron::execute(void) {
    systime = count * period; // time in seconds
    for (int i = 0; i < steps; ++i){
        solve(period / steps, y); // period in s
	output(0) = V * 1e-3; // convert to V
    } 
    count++;
}

void Alison_Neuron::update(DefaultGUIModel::update_flags_t flag) {
    switch (flag) {
        case INIT:
            setParameter("V0 (mV)", QString::number(V0)); 
            setParameter("cm (nF)", QString::number(cm)); 
            setParameter("G_Na_max (uS)", QString::number(G_Na_max)); 
            setParameter("E_Na (mV)", QString::number(E_Na)); 
            setParameter("G_K_max (uS)", QString::number(G_K_max)); 
            setParameter("G_Ks_max (uS)", QString::number(G_Ks_max)); 
            setParameter("E_K (mV)", QString::number(E_K)); 
            setParameter("G_L (uS)", QString::number(G_L)); 
            setParameter("E_L (mV)", QString::number(E_L)); 
            setParameter("G_A_max (uS)", QString::number(G_A_max)); 
            setParameter("Iapp (nA)", QString::number(Iapp)); 
            setParameter("Rate (Hz)", rate);
            setParameter("taua", QString::number(taua));
	        setParameter("taub1", QString::number(taub1));
            setState("m", m);
            setState("h", h);
            setState("n", n);
            setState("ns", ns);
            setState("b", b);
            setState("a", a);
            setState("IA",IA);
            setState("Time (s)", systime);
            break;

        case MODIFY:
            V0 = getParameter("V0 (mV)").toDouble();
            cm = getParameter("cm (nF)").toDouble();
            G_Na_max = getParameter("G_Na_max (uS)").toDouble();
            E_Na = getParameter("E_Na (mV)").toDouble();
            G_K_max = getParameter("G_K_max (uS)").toDouble();
            G_Ks_max = getParameter("G_Ks_max (uS)").toDouble();
            E_K = getParameter("E_K (mV)").toDouble();
            G_L = getParameter("G_L (uS)").toDouble();
            E_L = getParameter("E_L (mV)").toDouble();
            G_A_max = getParameter("G_A_max (uS)").toDouble();
            Iapp = getParameter("Iapp (nA)").toDouble(); 
            rate = getParameter("Rate (Hz)").toDouble();
            taua = getParameter("taua").toDouble();
            taub1 = getParameter("taub1").toDouble();
            steps = static_cast<int> (ceil(period * rate));

            break;

        case PERIOD:
            period = RT::System::getInstance()->getPeriod() * 1e-6; // time in seconds
            steps = static_cast<int> (ceil(period * rate));
            break;

        default:
            break;
    }
}

void Alison_Neuron::initParameters() {
    V0 = -55.038; // mV
    G_Na_max = 0.5;
    G_K_max = 0.02;
    G_Ks_max = 0.0;
    G_L = 0.001;
    G_A_max = 0.1;
    E_Na = 50.0; // mV
    E_K = -90.0;
    E_L = -59.5;
    cm = 0.0187;
    Iapp = 0.0; // 1 Hz spiking
    rate = 400;
    taua = 2.0;
    taub1 = 200.0;
    V = V0;
    m = 0.00001;
    h = 0.8522;
    n = 0.000208;
    ns = 0.0;
    a = 0.0;
    b = 0.0;
    count = 0;
    systime = 0;
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    steps = static_cast<int> (ceil(period * rate)); // calculate how many integrations to perform per execution step
}

void Alison_Neuron::solve(double dt, double *y) {
    double dydt[7];
    derivs(y, dydt);
    for (size_t i = 0; i < 7; ++i){
        y[i] += dt * dydt[i];
    }
}

void Alison_Neuron::derivs(double *y, double *dydt) {

    if(G_A_max == 0){
	IA = input(0) * (V - E_K); 
    }
    else{
	IA = G_A * (V - E_K)/cm;
    }	
    dV = -(G_Na*(V-E_Na) + G_K*(V-E_K) + G_Ks*(V-E_K) + G_L*(V-E_L))/cm + Iapp -IA + input(1);
    dm = (m_inf(V) - m) / tau_m(V);
    dh = (h_inf(V) - h) / tau_h(V);
    dn = (n_inf(V) - n) / tau_n(V);
    dns = (n_inf(V)-ns) / tau_ns(V);
    db = (b_inf(V) - b) / tau_b(V,taub1);
    da = (a_inf(V) - a) / taua;
}
