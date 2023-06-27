#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cpt.h"

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

static CPT<3, SAMPLE_RATE, FREQ> cpt;

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
	u[1] = in[1];
	u[2] = in[2];

	decltype(cpt)::power_vector i;

	i[0] = in[2];
	i[1] = in[3];
	i[2] = in[4];

	auto r = cpt.feed(u, i);
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
/** Test routine so this file can be used outside PSIM. **/
int main(int argc, char **argv) {
	static CPT<3, SAMPLE_RATE, FREQ> cpt;

	decltype(cpt)::power_vector u;
	decltype(cpt)::power_vector i;

	u.fill(2);
	i.fill(3);

	const auto r = cpt.feed(u, i);

	return 0;
}
#endif
