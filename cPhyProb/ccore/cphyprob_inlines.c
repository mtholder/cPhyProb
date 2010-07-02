/* We include this file in cphyprobs.h so that inlined functions will be inlined,
	but we also include this file as a compilation unit (for c compilers that
	don't support inline).

	This means that we have to use an odd system of sentinels to avoid:
		- infinite inclusion loops,
		- undefined symbols, and
		- multiply defined symbols
*/
#if !defined(C_PHY_PROB_INLINE)
#define C_PHY_PROB_INLINE


#if defined(C_PHY_PROB_H)
#	if defined(NO_INLINE) && NO_INLINE
#		define DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT 0
#	else
#		define DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT 1
#	endif
#else
#	include "cphyprob.h"
#	if defined(NO_INLINE) && NO_INLINE
#		define DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT 1
#	else
#		define DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT 0
#	endif
#endif

#if defined(NO_INLINE) && NO_INLINE
#	define INLINE
#else
#	define INLINE inline
#endif

#if defined(DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT) && DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT
/*
#	if defined(USE_BEAGLE_LIB) && USE_BEAGLE_LIB
#		include "libhmsbeagle/beagle.h"
		INLINE void BEAGLE_ERROR_MSG(int rc, const char * funcName) {
				if (rc != BEAGLE_SUCCESS) {
					CPHYPROB_DEBUG_PRINTF2("Beagle error code %d in function \"%s\"", rc, funcName);
				}
			}
#	endif
*/

	/* End of odd sentinel */
#	undef DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT
#	endif //defined(DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT) && DO_INCLUDED_C_PHY_PROB_INLINE_CONTENT
#endif //(C_PHY_PROB_INLINE)
