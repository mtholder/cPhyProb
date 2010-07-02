#include "state_set_lookup.h"


void state_set_lookup_dtor(StateSetLookupStruct* ssl_struct) {
	CPHYPROB_DEBUG_PRINTF("In state_set_lookup_dtor\n");
	if (ssl_struct == 0L)
		return;
	if (ssl_struct->state_lookup) {
		if (ssl_struct->state_lookup[0])
			free(ssl_struct->state_lookup[0]);
		free(ssl_struct->state_lookup);
	}
	PyObject_Del(ssl_struct);
	CPHYPROB_DEBUG_PRINTF("Leaving state_set_lookup_dtor\n");
}



