#include <stdlib.h>
#include <string.h>

#include "../cpt.h"

#ifdef _WIN32
#define ADDAPI __declspec(dllexport)
#define ADDCALL __cdecl
#else
#define ADDAPI
#define ADDCALL
#endif

#define SAMPLE_RATE	12000
#define FREQ        60

#ifdef _WIN32

#ifdef __cplusplus
extern "C"
{
#endif

// PLACE GLOBAL VARIABLES OR USER FUNCTIONS HERE...

static CPT<1, SAMPLE_RATE, FREQ> cpt;

/////////////////////////////////////////////////////////////////////
// FUNCTION: SimulationStep
//   This function runs at every time step.
//double t: (read only) time
//double delt: (read only) time step as in Simulation control
//double *in: (read only) zero based array of input values. in[0] is the first node, in[1] second input...
//double *out: (write only) zero based array of output values. out[0] is the first node, out[1] second output...
//int *pnError: (write only)  assign  *pnError = 1;  if there is an error and set the error message in szErrorMsg
//    strcpy(szErrorMsg, "Error message here..."); 
// DO NOT CHANGE THE NAME OR PARAMETERS OF THIS FUNCTION
ADDAPI void ADDCALL SimulationStep(
		double t, double delt, double *in, double *out,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr) {
// ENTER YOUR CODE HERE...

	decltype(cpt)::power_vector u;

	u[0] = in[0];

	decltype(cpt)::power_vector i;

	i[0] = in[1];

	const auto r = cpt.feed(u, i);

	out[0] = r.u[0];

	out[1] = r.i[0];

	out[2] = r.p[0];

	out[3] = r.P[0];

	out[4] = r.w[0];

	out[5] = r.W[0];

	out[6] = r.U[0];

	out[7] = r.ia[0];

	out[8] = r.ir[0];

	out[9] = r.iv[0];

	out[10] = r.iba[0];

	out[11] = r.ibr[0];

	out[12] = r.iua[0];

	out[13] = r.iur[0];

	out[14] = r.t.P;

	out[15] = r.t.Q;

	out[16] = r.t.V;

	out[17] = r.t.A;

	out[18] = r.t.rf;

	out[19] = r.t.uf;

	out[20] = r.t.df;

	out[21] = r.t.pf;
}


/////////////////////////////////////////////////////////////////////
// FUNCTION: SimulationBegin
//   Initialization function. This function runs once at the beginning of simulation
//   For parameter sweep or AC sweep simulation, this function runs at the beginning of each simulation cycle.
//   Use this function to initialize static or global variables.
//const char *szId: (read only) Name of the C-block 
//int nInputCount: (read only) Number of input nodes
//int nOutputCount: (read only) Number of output nodes
//int nParameterCount: (read only) Number of parameters is always zero for C-Blocks.  Ignore nParameterCount and pszParameters
//int *pnError: (write only)  assign  *pnError = 1;  if there is an error and set the error message in szErrorMsg
//    strcpy(szErrorMsg, "Error message here..."); 
// DO NOT CHANGE THE NAME OR PARAMETERS OF THIS FUNCTION
ADDAPI void ADDCALL SimulationBegin(
		const char *szId, int nInputCount, int nOutputCount,
		 int nParameterCount, const char ** pszParameters,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr) {
// ENTER INITIALIZATION CODE HERE...




}


/////////////////////////////////////////////////////////////////////
// FUNCTION: SimulationEnd
//   Termination function. This function runs once at the end of simulation
//   For parameter sweep or AC sweep simulation, this function runs at the end of each simulation cycle.
//   Use this function to de-allocate any allocated memory or to save the result of simulation in an alternate file.
// Ignore all parameters for C-block 
// DO NOT CHANGE THE NAME OR PARAMETERS OF THIS FUNCTION
ADDAPI void ADDCALL SimulationEnd(const char *szId, void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr) {




}

#ifdef __cplusplus
} // __cplusplus defined.
#endif

#else

#include <fstream>

/** Test routine so this file can be used outside PSIM. **/
int main(int argc, char **argv) {
	static CPT<1, SAMPLE_RATE, FREQ> cpt;

	decltype(cpt)::power_vector u;
	decltype(cpt)::power_vector i;

	std::ofstream file ("/tmp/.marcos-private/saida.csv", std::ofstream::out);

	file << "Va" << "\t" << "Vb" << "\t" << "Vc" << "\t";
	file << "Ia" << "\t" << "Ib" << "\t" << "Ic" << "\t";
	file << "pa" << "\t" << "pb" << "\t" << "pc" << "\t";
	file << "Pa" << "\t" << "Pb" << "\t" << "Pc" << "\t";
	file << "wa" << "\t" << "wb" << "\t" << "wc" << "\t";
	file << "Wa" << "\t" << "Wb" << "\t" << "Wc" << "\t";
	file << "Ua" << "\t" << "Ub" << "\t" << "Uc" << "\t";
	file << "iaa" << "\t" << "iab" << "\t" << "iac" << "\t";
	file << "ira" << "\t" << "irb" << "\t" << "irc" << "\t";
	file << "iva" << "\t" << "ivb" << "\t" << "ivc";
	file << std::endl;

	const auto angle = -M_PI/2;

	for (unsigned int o=0; o<8*SAMPLE_RATE/FREQ; o++) {
		u[0] =   1 * M_SQRT2 * sin(2 * M_PI * o / (SAMPLE_RATE/FREQ) );
		u[1] =   0 * M_SQRT2 * sin(2 * M_PI * o / (SAMPLE_RATE/FREQ) - 2*M_PI/3);
		u[2] =   0 * M_SQRT2 * sin(2 * M_PI * o / (SAMPLE_RATE/FREQ) + 2*M_PI/3);

		i[0] =   1 * M_SQRT2 * sin(2 * M_PI * o / (SAMPLE_RATE/FREQ) + angle);
		i[1] =   0 * M_SQRT2 * sin(2 * M_PI * o / (SAMPLE_RATE/FREQ) - 2*M_PI/3 + angle);
		i[2] =   0 * M_SQRT2 * sin(2 * M_PI * o / (SAMPLE_RATE/FREQ) + 2*M_PI/3 + angle);

		const auto r = cpt.feed(u, i);

		file << r.u[0] << "\t" << r.u[1] << "\t" << r.u[2] << "\t";
		file << r.i[0] << "\t" << r.i[1] << "\t" << r.i[2] << "\t";
		file << r.p[0] << "\t" << r.p[1] << "\t" << r.p[2] << "\t";
		file << r.P[0] << "\t" << r.P[1] << "\t" << r.P[2] << "\t";
		file << r.w[0] << "\t" << r.w[1] << "\t" << r.w[2] << "\t";
		file << r.W[0] << "\t" << r.W[1] << "\t" << r.W[2] << "\t";
		file << r.U[0] << "\t" << r.U[1] << "\t" << r.U[2] << "\t";
		file << r.ia[0] << "\t" << r.ia[1] << "\t" << r.ia[2] << "\t";
		file << r.ir[0] << "\t" << r.ir[1] << "\t" << r.ir[2] << "\t";
		file << r.iv[0] << "\t" << r.iv[1] << "\t" << r.iv[2];
		file << std::endl;
	}

	return 0;
}
#endif
