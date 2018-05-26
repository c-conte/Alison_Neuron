#include <default_gui_model.h>


class Alison_Neuron : public DefaultGUIModel
{
    
    public:
    
        Alison_Neuron(void);
        virtual~Alison_Neuron(void);
    
        void execute(void);
    
    protected:
    
        void update(DefaultGUIModel::update_flags_t);
    
    private:
    
        void derivs(double *, double *);
        void solve(double, double *);
        void initParameters();
    
        double y[7];
        double period;
        int steps;
    
        double V0;
        double cm;
        double G_Na_max;
        double E_Na;
        double G_K_max;
        double E_K;
        double G_L;
        double E_L;
        double G_A_max;
        double Iapp;
        double rate;
        double taua;
	double taub;
	double taun0;
	double taun1;
	double tauh0;
	double tauh1;
        double IA;
        double systime;
        long long count;
};
