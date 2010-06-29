/**
 *
 * Extension adapted from Alex Martelli's example code in the Python Cookbook 17.1
 *
 * ASRV code from the PAML package (ref below):
 *
 *
 * Yang, Z.  2007.  PAML 4: a program package for phylogenetic analysis by maximum
 *	likelihood.  Molecular Biology and Evolution 24: 1586-1591
 *	(http://abacus.gene.ucl.ac.uk/software/paml.html)
 *
 *
 *  Eigensystem code from MrBayes 3.2
 *  by Fredrik Ronquist, John P. Huelsenbeck, and Paul van der Mark
 *  updated from the CVS Repository: mrbayes from the root:
 *   :pserver:anonymous@mrbayes.cvs.sourceforge.net:/cvsroot/mrbayes
 *  on 2007-Oct-05
 *
 *  Copyright 2002-2007
 *  	John P. Huelsenbeck
 *  	Fredrik Ronquist
 *  	Paul van der Mark
 *
 *	Some of that code was written (or translated from FORTRAN) by David Swofford
 *
 * Other code Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
 * mtholder@gmail.com
 * (see bottom of file)
 *
 *
 */
#include "dsct_model.h"
#define ALL_IN_ONE

#if defined(ALL_IN_ONE)
#	include "dsct_model.c"
#else
	extern const int UNDERFLOW_ERROR_CODE;
#endif
static PyObject * CLAUnderflowError;


/* python C API utility funcs */
static DSCTModelObj *accessor_args_to_DSCTModel(PyObject *self, PyObject *args);
static DSCTModelObj *accessor_args_to_DSCTModel_i(PyObject *self, PyObject *args, int *index);
static DSCTModelObj *accessor_args_to_DSCTModel_i_i(PyObject *self, PyObject *args, int *r, int *c);
static DSCTModelObj *accessor_args_to_DSCTModel_i_i_d(PyObject *self, PyObject *args, int *r, int *c, double *d);
static PMatArrayObj *accessor_args_to_PMat(PyObject *self, PyObject *args);
static PyObject * doubleArrayToList(const double *arr, unsigned len);
static PyObject * doubleMatToList(const double **arr, unsigned n_rows, unsigned n_cols);
static PyObject * double3DMatToList(const double ***arr, unsigned n_mats, unsigned n_rows, unsigned n_cols);
static PyObject * listToDoubleArray(PyObject *list_obj, double *arr, unsigned n, int demandExactLen);
static PyObject * listToDoubleMatrix(PyObject *list_obj, double **arr, unsigned n_rows, unsigned n_cols, int demandExactLen);
static PyObject* none();
static int pyIndToCInd(int py_ind, unsigned dim);
/*Exposed function prototypes*/
static PyObject* csslookup_ctor(PyObject *self, PyObject *args);
static PyObject* cleaf_data_ctor(PyObject *self, PyObject *args);
static PyObject* cla_ctor(PyObject *self, PyObject *args);
static PyObject* full_la_ctor(PyObject *self, PyObject *args);
static PyObject* cdsctm_ctor(PyObject *self, PyObject *args);
static PyObject* cdsctm_get_dimen(PyObject *self, PyObject *args);
static PyObject* cdsctm_get_q_mat(PyObject *self, PyObject *args);
static PyObject* cdsctm_get_q_mat_row(PyObject *self, PyObject *args);
static PyObject* cdsctm_get_q_mat_cell(PyObject *self, PyObject *args);
static PyObject* cdsctm_set_q_mat(PyObject *self, PyObject *args);
static PyObject* cdsctm_set_q_mat_row(PyObject *self, PyObject *args);
static PyObject* cdsctm_set_q_mat_cell(PyObject *self, PyObject *args);
static PyObject* cpmat_array_ctor(PyObject *self, PyObject *args);
static PyObject* cpmat_array_get_n_matrices(PyObject *self, PyObject *args);
static PyObject* cpmat_array_get_n_states(PyObject *self, PyObject *args);
static PyObject* cpmat_array_get_mat_list(PyObject *self, PyObject *args);
static PyObject* cpmat_array_get_mat(PyObject *self, PyObject *args);
static PyObject* casrvo_ctor(PyObject *self, PyObject *args);
static PyObject* casrvo_get_n_cat(PyObject *self, PyObject *args);
static PyObject* casrvo_get_rates(PyObject *self, PyObject *args);
static PyObject* casrvo_get_shape(PyObject *self, PyObject *args);
static PyObject* casrvo_set_shape(PyObject *self, PyObject *args);


/**
 * (internal) converts an python index on [-dim,dim) to the positive version
 * of the index.
 * Returns a non-negative index for success or -1.
 * if -1 is returned then the function will have generated a
 * IndexError("list index out of range") already.
 *
 */
static int pyIndToCInd(int py_ind, unsigned dim) {
	PRINTF2("In pyIndToCInd(%d,%u)\n", py_ind, dim);
	if (py_ind < 0)
		py_ind += (int)dim;
	if (py_ind < (int)dim)
		return py_ind;
	PyErr_SetString(PyExc_IndexError, "list index out of range");
	return -1;
}
/**
 * (internal) Converts a double array to a python list of floats
 *
 *	\assert(arr != 0L || len == 0)
 */
static PyObject * doubleArrayToList(const double *arr, unsigned len) {
	PyObject *lp;
	PyObject *el_obj;
	unsigned i;
	assert(arr != 0L || len == 0);
	lp = PyList_New(len);
	if (lp == 0L)
		return 0L;
	for (i = 0; i < len; ++i) {
		el_obj = PyFloat_FromDouble(arr[i]);
		if (el_obj == 0L) {
			Py_DECREF(lp);
			return 0L;
		}
		PyList_SetItem(lp, i, el_obj);
	}
	return lp;
}
/**
 * (internal)  Converts a rectangular 2D matrix of floats to a python list of
 * 	list of floats.
 *
 *	\assert(arr != 0L || len == 0)
 */
static PyObject * doubleMatToList(const double **arr, unsigned n_rows, unsigned n_cols) {
	PyObject *lp;
	PyObject *el_obj;
	unsigned i;
	assert(arr != 0L || n_rows == 0);
	lp = PyList_New(n_rows);
	if (lp == 0L)
		return 0L;
	for (i = 0; i < n_rows; ++i) {
		el_obj = doubleArrayToList(arr[i], n_cols);
		if (el_obj == 0L) {
			Py_DECREF(lp);
			return 0L;
		}
		PyList_SetItem(lp, i, el_obj);
	}
	return lp;
}
/**
 * (internal)  Converts a rectangular 3D matrix of floats to a python list of
 * 	list of floats.
 *
 *	\assert(arr != 0L || len == 0)
 */
static PyObject * double3DMatToList(const double ***arr, unsigned n_mats, unsigned n_rows, unsigned n_cols) {
	PyObject *lp;
	PyObject *el_obj;
	unsigned i;
	assert(arr != 0L || n_mats == 0);
	lp = PyList_New(n_mats);
	if (lp == 0L)
		return 0L;
	for (i = 0; i < n_mats; ++i) {
		el_obj = doubleMatToList(arr[i], n_rows, n_cols);
		if (el_obj == 0L) {
			Py_DECREF(lp);
			return 0L;
		}
		PyList_SetItem(lp, i, el_obj);
	}
	return lp;
}
static PyObject * listToDoubleArray(PyObject *list_obj, double *arr, unsigned n, int demandExactLen) {
	PyObject *item, *f_item;
	unsigned pylist_len = (unsigned) PyList_Size(list_obj);
	unsigned i;
	if (pylist_len != n && ((demandExactLen != 0) || (pylist_len < n))) {
		PyErr_SetString(PyExc_IndexError, "list index out of range");
		return 0L;
	}
	for (i = 0; i < n; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		f_item = PyNumber_Float(item);
		if (f_item == 0L) {
			Py_DECREF(item);
			return 0L;
		}
		arr[i] = PyFloat_AsDouble(item);
		Py_DECREF(item);
		Py_DECREF(f_item);
	}
	return none();
}
static PyObject * listToDoubleMatrix(PyObject *list_obj, double **arr, unsigned n_rows, unsigned n_cols, int demandExactLen) {
	PyObject *item, *r_item;
	unsigned pylist_len = (unsigned) PyList_Size(list_obj);
	unsigned i;
	if (pylist_len != n_rows && ((demandExactLen != 0) || (pylist_len < n_rows))) {
		PyErr_SetString(PyExc_IndexError, "list index out of range");
		return 0L;
	}
	for (i = 0; i < n_rows; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		if (!PyList_Check(item)) {
			PyErr_SetString(PyExc_TypeError, "list of list of floats expected");
			Py_DECREF(item);
			return 0L;
		}
		r_item = listToDoubleArray(item, arr[i], n_cols, demandExactLen);
		Py_DECREF(item);
		if (0L == r_item)
			return 0L;
	}
	return none();
}


/* A minimal Python type-object */
statichere PyTypeObject state_set_lookup_type = {
	PyObject_HEAD_INIT(0)	  /* initialize to 0 to ensure Win32 portability  */
	0,						  /* ob_size */
	"state_set_lookup",			 /* tp_name */
	sizeof(StateSetLookupStruct),		/* tp_basicsize */
	0,						  /* tp_itemsize */
	/* methods */
	(destructor)sslookup_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
};
statichere PyTypeObject leaf_data_type = {
	PyObject_HEAD_INIT(0)	  /* initialize to 0 to ensure Win32 portability  */
	0,						  /* ob_size */
	"leaf_data",					/* tp_name */
	sizeof(LeafDataObj),		/* tp_basicsize */
	0,						  /* tp_itemsize */
	/* methods */
	(destructor)leaf_data_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
};
statichere PyTypeObject cla_type = {
	PyObject_HEAD_INIT(0)	  /* initialize to 0 to ensure Win32 portability  */
	0,						  /* ob_size */
	"cla",						/* tp_name */
	sizeof(CLAObj),		/* tp_basicsize */
	0,						  /* tp_itemsize */
	/* methods */
	(destructor)cla_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
};
statichere PyTypeObject full_la_type = {
	PyObject_HEAD_INIT(0)	  /* initialize to 0 to ensure Win32 portability  */
	0,						  /* ob_size */
	"full_la",						/* tp_name */
	sizeof(FullLAObj),		/* tp_basicsize */
	0,						  /* tp_itemsize */
	/* methods */
	(destructor)full_la_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
};
statichere PyTypeObject dsct_model_type = {
	PyObject_HEAD_INIT(0)	  /* initialize to 0 to ensure Win32 portability  */
	0,						  /* ob_size */
	"dsct_model",					/* tp_name */
	sizeof(DSCTModelObj),		/* tp_basicsize */
	0,						  /* tp_itemsize */
	/* methods */
	(destructor)cdsctm_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
};
statichere PyTypeObject pmat_array_type = {
	PyObject_HEAD_INIT(0)	  /* initialize to 0 to ensure Win32 portability  */
	0,						  /* ob_size */
	"pmat_array",					/* tp_name */
	sizeof(PMatArrayObj),		/* tp_basicsize */
	0,						  /* tp_itemsize */
	/* methods */
	(destructor)cpmat_array_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
};
statichere PyTypeObject asrv_type = {
	PyObject_HEAD_INIT(0)	  	/* initialize to 0 to ensure Win32 portability  */
	0,						  	/* ob_size */
	"asrv_obj",	/* tp_name */
	sizeof(ASRVObj),		/* tp_basicsize */
	0,						  	/* tp_itemsize */
	/* methods */
	(destructor)asrv_obj_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
};
static DSCTModelObj *accessor_args_to_DSCTModel(PyObject *self, PyObject *args) {
	PyObject *dsct_model_py_obj;
	if (!PyArg_ParseTuple(args, "O!", &dsct_model_type, &dsct_model_py_obj))
		return 0L;
	return (DSCTModelObj *)(dsct_model_py_obj);
}
static DSCTModelObj *accessor_args_to_DSCTModel_i(PyObject *self, PyObject *args, int *index) {
	PyObject *dsct_model_py_obj;
	if (!PyArg_ParseTuple(args, "O!i", &dsct_model_type, &dsct_model_py_obj, index))
		return 0L;
	return (DSCTModelObj *)(dsct_model_py_obj);
}
static DSCTModelObj *accessor_args_to_DSCTModel_i_i(PyObject *self, PyObject *args, int *r, int *c) {
	PyObject *dsct_model_py_obj;
	if (!PyArg_ParseTuple(args, "O!ii", &dsct_model_type, &dsct_model_py_obj, r, c))
		return 0L;
	return (DSCTModelObj *)(dsct_model_py_obj);
}
static DSCTModelObj *accessor_args_to_DSCTModel_i_i_d(PyObject *self, PyObject *args, int *r, int *c, double *d) {
	PyObject *dsct_model_py_obj;
	if (!PyArg_ParseTuple(args, "O!iid", &dsct_model_type, &dsct_model_py_obj, r, c, d))
		return 0L;
	return (DSCTModelObj *)(dsct_model_py_obj);
}
static PMatArrayObj *accessor_args_to_PMat(PyObject *self, PyObject *args) {
	PyObject *pmat_py_obj;
	if (!PyArg_ParseTuple(args, "O!", &pmat_array_type, &pmat_py_obj))
		return 0L;
	return (PMatArrayObj *)(pmat_py_obj);
}
static PyObject* none() {
	Py_INCREF(Py_None);
	return Py_None;
}
/* MODULE functions */
/* the factory function */
static PyObject* csslookup_ctor(PyObject *self, PyObject *args) {
	PyObject * list_obj = 0L;
	PyObject *item, *f_item, *s_item;
	StateSetLookupStruct * to_return = 0L;
	int n_states_declared = 0;
	unsigned n_states = 0;
	unsigned n_state_sets = 0;
	unsigned arr_len = 0;
	unsigned row_len = 0;
	unsigned i;
	int j;
	int curr_row_len, curr_el;
	int **p = 0L;
	int *curr_p;
	if (!PyArg_ParseTuple(args, "iO!", &n_states_declared, &PyList_Type, &list_obj))
		return 0;
	if (n_states_declared < 2) {
		PyErr_SetString(PyExc_ValueError, "The number of states must be greater than 1.");
		return 0;
	}
	n_states = (unsigned) n_states_declared;
	n_state_sets = (unsigned) PyList_Size(list_obj);
	if (n_state_sets < 2) {
		PyErr_SetString(PyExc_IndexError, "state set lookup must have at least two states.");
		return 0L;
	}
	for (i = 0; i < n_state_sets; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		if (!PyList_Check(item)) {
			PyErr_SetString(PyExc_TypeError, "list of list of ints expected");
			Py_DECREF(item);
			return 0L;
		}
		arr_len += 1; //add one for the # of states element
		row_len = (unsigned) PyList_Size(item);
		if (row_len < 1) {
			PyErr_SetString(PyExc_IndexError, "each state set must have at least one states.");
			return 0L;
		}
		if (row_len > n_states_declared) {
			PyErr_SetString(PyExc_IndexError, "The state set cannot contain more states than the n_states.");
			return 0L;
		}
		arr_len += row_len;
	}
	/*allocate a ragged two-D array with the memory contiguous.*/
	p = (int**)malloc(n_state_sets*sizeof(int*));
	if (p == 0L) {
		PyErr_NoMemory();
		return 0L;
	}
	p[0] = (int*)malloc(arr_len*sizeof(int));
	if (p[0] == 0L) {
		PyErr_NoMemory();
		return 0L;
	}
	arr_len = 0;
	curr_p = p[0];
	for (i = 0; i < n_state_sets; ++i) {
		p[i] = curr_p;
		item = PyList_GetItem(list_obj, i);
		assert(item != 0L);
		Py_INCREF(item);
		if (!PyList_Check(item)) {
			PyErr_SetString(PyExc_TypeError, "list of list of ints expected");
			Py_DECREF(item);
			return 0L;
		}
		arr_len += 1; //add one for the # of states element
		curr_row_len = (int) PyList_Size(item);
		assert(curr_row_len > 0);
		*curr_p++ = curr_row_len;
		for (j = 0; j < curr_row_len; ++j) {
			f_item = PyList_GetItem(item, j);
			if (f_item == 0L)
				goto errorExit;
			Py_INCREF(f_item);
			s_item = PyNumber_Int(f_item);
			if (s_item == 0L) {
				Py_DECREF(f_item);
				goto errorExit;
			}
			curr_el = PyInt_AsLong(f_item);
			Py_DECREF(f_item);
			Py_DECREF(s_item);
			*curr_p++ = curr_el;
		}
	}
	/*donate the pointer to the object we are creating.*/
	to_return = sslookup_new(n_states, n_state_sets, p);
	if (to_return == 0L)
		goto errorExit;
	return (PyObject*) to_return;
	errorExit:
		PRINTF("In csslookup_ctor errorExit\n");
		free(p[0]);
		free(p);
		return 0L;
}
static PyObject* cleaf_data_ctor(PyObject *self, PyObject *args) {
	PRINTF("In cleaf_data_ctor\n");
	PyObject * list_obj = 0L;
	PyObject * sslookup_obj = 0L;
	PyObject * item, *f_item;
	StateSetLookupStruct * sslookup_ptr;
	LeafDataObj * to_return;
	unsigned n_sites = 0;
	int n_cats_int = 0;
	unsigned n_cats;
	int * curr_p;
	int curr_el, i;
	if (!PyArg_ParseTuple(args, "O!O!i", &PyList_Type, &list_obj, &state_set_lookup_type, &sslookup_obj, &n_cats_int)){
		return 0L;
	}
	sslookup_ptr = (StateSetLookupStruct *) sslookup_obj;	
	if (n_cats_int < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of model categories must be greater than 0.");
		return 0;
	}
	n_cats = (unsigned) n_cats_int;
	n_sites = (unsigned) PyList_Size(list_obj);
	if (n_sites < 1) {
		PyErr_SetString(PyExc_IndexError, "n_sites must be greater than 0");
		return 0L;
	}
	to_return = leaf_data_new(n_sites, n_cats, sslookup_ptr);
	if (to_return == 0L)
		return 0L;
	curr_p = to_return->ssind;
	for (i = 0; i < n_sites; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L)
			goto errorExit;
		Py_INCREF(item);
		f_item = PyNumber_Int(item);
		if (f_item == 0L) {
			Py_DECREF(item);
			goto errorExit;
		}
		curr_el = PyInt_AsLong(item);
		Py_DECREF(item);
		Py_DECREF(f_item);
		*curr_p++ = curr_el;
	}
	return (PyObject*) to_return;
	errorExit:
		PRINTF("In cleaf_data_ctor errorExit\n");
		leaf_data_dtor(to_return);
		return 0L;
}
static PyObject* cla_ctor(PyObject *self, PyObject *args) {
	int n_sites_int = 0;
	unsigned n_sites = 0;
	int n_states_int = 0;
	unsigned n_states = 0;
	int n_cats_int = 0;
	unsigned n_cats;
	if (!PyArg_ParseTuple(args, "iii", &n_sites_int, &n_states_int, &n_cats_int))
		return 0;
	if (n_cats_int < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of model categories must be greater than 0.");
		return 0;
	}
	if (n_sites_int < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of sites must be greater than 0.");
		return 0;
	}
	if (n_states_int < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of states must be greater than 0.");
		return 0;
	}
	n_cats = (unsigned) n_cats_int;
	n_sites = (unsigned) n_sites_int;
	n_states = (unsigned) n_states_int;
	return (PyObject*) cla_new(n_sites, n_states, n_cats);
}
static PyObject* full_la_set_freqs(PyObject *self, PyObject *args) {
	PyObject * full_la_py_obj;
	PyObject * categ_freq_py_obj;
	PyObject * state_freq_py_obj;
	FullLAObj * full_la;
	unsigned n_categ, n_states;
	PyObject *item, *f_item, *to_return;
	unsigned pylist_len;
	unsigned i;
	
	if (!PyArg_ParseTuple(args, "O!O!O!", &full_la_type, &full_la_py_obj, &PyList_Type, &categ_freq_py_obj, &PyList_Type, &state_freq_py_obj))
		return 0L;
	full_la = (FullLAObj *)(full_la_py_obj);
	if (full_la->cso.n_categ_arr) {
		PyErr_SetString(PyExc_ValueError, "full_la_set_freqs should not be called when using paritioned models");
		return 0L;
	}
	n_categ = full_la->cso.n_categ_or_subs;
	n_states = full_la->n_states;
	
	pylist_len  = (unsigned) PyList_Size(categ_freq_py_obj);
	if (pylist_len != n_categ) {
		PyErr_SetString(PyExc_IndexError, "Expecting a category frequency for every category");
		return 0L;
	}
	for (i = 0; i < n_categ; ++i) {
		item = PyList_GetItem(categ_freq_py_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		f_item = PyNumber_Float(item);
		if (f_item == 0L) {
			Py_DECREF(item);
			return 0L;
		}
		full_la->state_categ_freqs[i][n_states] = PyFloat_AsDouble(item);
		Py_DECREF(item);
		Py_DECREF(f_item);
	}
	to_return = listToDoubleMatrix(state_freq_py_obj, full_la->state_categ_freqs, n_categ, n_states, 1);
	return to_return;
}
static PyObject* full_la_ctor(PyObject *self, PyObject *args) {
	int n_sites_int = 0;
	unsigned n_sites = 0;
	int n_states_int = 0;
	unsigned n_states = 0;
	int n_cats_int = 0;
	unsigned n_cats;
	if (!PyArg_ParseTuple(args, "iii", &n_sites_int, &n_states_int, &n_cats_int))
		return 0;
	if (n_cats_int < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of model categories must be greater than 0.");
		return 0;
	}
	if (n_sites_int < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of sites must be greater than 0.");
		return 0;
	}
	if (n_states_int < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of states must be greater than 0.");
		return 0;
	}
	n_cats = (unsigned) n_cats_int;
	n_sites = (unsigned) n_sites_int;
	n_states = (unsigned) n_states_int;
	return (PyObject*) full_la_new(n_sites, n_states, n_cats);
}
static PyObject* cdsctm_ctor(PyObject *self, PyObject *args) {
	int dimen;
	unsigned dim;
	if (!PyArg_ParseTuple(args, "i", &dimen))
		return 0;
	if (dimen < 2) {
		PyErr_SetString(PyExc_ValueError, "The dimensions parameter must be greater than 1");
		return 0;
	}
	dim = (unsigned) dimen;
	return (PyObject*)dsct_model_new(dim);
}
/* dimensions (n_states) accessor */
static PyObject* cdsctm_get_dimen(PyObject *self, PyObject *args) {
	DSCTModelObj *dsct_model_obj = accessor_args_to_DSCTModel(self, args);
	if (dsct_model_obj == 0L)
		return 0L;
	return PyInt_FromLong((long)(dsct_model_obj->dim));
}
/* Q matrix accessor */
static PyObject* cdsctm_get_q_mat(PyObject *self, PyObject *args) {
	DSCTModelObj *dsct_model_obj = accessor_args_to_DSCTModel(self, args);
	if (dsct_model_obj == 0L)
		return 0L;
	return doubleMatToList((const double **)dsct_model_obj->q_mat, dsct_model_obj->dim, dsct_model_obj->dim);
}
/* Q matrix row accessor */
static PyObject* cdsctm_get_q_mat_row(PyObject *self, PyObject *args) {
	int row_ind;
	unsigned dim;
	DSCTModelObj *dsct_model_obj = accessor_args_to_DSCTModel_i(self, args, &row_ind);
	if (dsct_model_obj == 0L)
		return 0L;
	dim = dsct_model_obj->dim;
	row_ind  = pyIndToCInd(row_ind, dim);
	if (row_ind < 0) {
		return 0L;
	}
	assert(dsct_model_obj->q_mat);
	assert(dsct_model_obj->q_mat[row_ind]);
	return doubleArrayToList(dsct_model_obj->q_mat[row_ind], dim);
}
/* Q matrix cell accessor */
static PyObject* cdsctm_get_q_mat_cell(PyObject *self, PyObject *args) {
	int row_ind, col_ind;
	unsigned dim;
	DSCTModelObj *dsct_model_obj = accessor_args_to_DSCTModel_i_i(self, args, &row_ind, &col_ind);
	if (dsct_model_obj == 0L)
		return 0L;
	dim = dsct_model_obj->dim;
	row_ind  = pyIndToCInd(row_ind, dim);
	if (row_ind < 0) {
		return 0L;
	}
	col_ind  = pyIndToCInd(col_ind, dim);
	if (col_ind < 0) {
		return 0L;
	}
	assert(dsct_model_obj->q_mat);
	assert(dsct_model_obj->q_mat[row_ind]);
	assert(dsct_model_obj->q_mat[row_ind][col_ind]);
	return PyFloat_FromDouble(dsct_model_obj->q_mat[row_ind][col_ind]);
}
/* Q matrix setter */
static PyObject* cdsctm_set_q_mat(PyObject *self, PyObject *args) {
	unsigned dim;
	PyObject * list_obj = 0L;
	PyObject * dsct_model_py_obj;
	PyObject * to_return;
	DSCTModelObj *dsct_model_obj;
	if (!PyArg_ParseTuple(args, "O!O!", &dsct_model_type, &dsct_model_py_obj, &PyList_Type, &list_obj))
		return 0L;
	dsct_model_obj = (DSCTModelObj *)(dsct_model_py_obj);
	dim = dsct_model_obj->dim;
	to_return = listToDoubleMatrix(list_obj, dsct_model_obj->q_mat, dim, dim, 1);
	dsct_model_obj->eigen_calc_dirty = 1;
	return to_return;
}
/* Q matrix row setter */
static PyObject* cdsctm_set_q_mat_row(PyObject *self, PyObject *args) {
	int row_ind;
	unsigned dim;
	PyObject * to_return;
	PyObject * list_obj = 0L;
	PyObject *dsct_model_py_obj;
	DSCTModelObj *dsct_model_obj;
	if (!PyArg_ParseTuple(args, "O!iO!", &dsct_model_type, &dsct_model_py_obj, &row_ind, &PyList_Type, &list_obj))
		return 0L;
	dsct_model_obj = (DSCTModelObj *)(dsct_model_py_obj);
	dim = dsct_model_obj->dim;
	row_ind  = pyIndToCInd(row_ind, dim);
	if (row_ind < 0)
		return 0L;
	to_return = listToDoubleArray(list_obj, dsct_model_obj->q_mat[row_ind], dim, 1);
	dsct_model_obj->eigen_calc_dirty = 1;
	return to_return;
}
/* Q matrix cell setter */
static PyObject* cdsctm_set_q_mat_cell(PyObject *self, PyObject *args) {
	int row_ind, col_ind;
	double val;
	unsigned dim;
	DSCTModelObj *dsct_model_obj = accessor_args_to_DSCTModel_i_i_d(self, args, &row_ind, &col_ind, &val);
	if (dsct_model_obj == 0L)
		return 0L;
	dim = dsct_model_obj->dim;
	row_ind  = pyIndToCInd(row_ind, dim);
	if (row_ind < 0) {
		return 0L;
	}
	col_ind  = pyIndToCInd(col_ind, dim);
	if (col_ind < 0) {
		return 0L;
	}
	assert(dsct_model_obj->q_mat);
	assert(dsct_model_obj->q_mat[row_ind]);
	assert(dsct_model_obj->q_mat[row_ind][col_ind]);
	dsct_model_obj->q_mat[row_ind][col_ind] = val;
	dsct_model_obj->eigen_calc_dirty = 1;
	return none();
}
static PyObject* cpmat_array_ctor(PyObject *self, PyObject *args){
	int n_mat, n_states;
	unsigned un_mat, un_states;
	if (!PyArg_ParseTuple(args, "ii", &n_mat, &n_states))
		return 0;
	if (n_mat < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of matrices must be greater than 0");
		return 0;
	}
	if (n_states < 2) {
		PyErr_SetString(PyExc_ValueError, "The number of states must be greater than 1");
		return 0;
	}
	un_states = (unsigned) n_states;
	un_mat = (unsigned) n_mat;
	return (PyObject*)cpmat_array_new(un_mat, un_states);
}
/* `n_mat` accessor */
static PyObject* cpmat_array_get_n_matrices(PyObject *self, PyObject *args) {
	PMatArrayObj *pmat_obj = accessor_args_to_PMat(self, args);
	if (pmat_obj == 0L)
		return 0L;
	return PyInt_FromLong((long)(pmat_obj->n_mat));
}
/* `n_states` accessor */
static PyObject* cpmat_array_get_n_states(PyObject *self, PyObject *args) {
	PMatArrayObj *pmat_obj = accessor_args_to_PMat(self, args);
	if (pmat_obj == 0L)
		return 0L;
	return PyInt_FromLong((long)(pmat_obj->n_states));
}
static PyObject* cpmat_array_get_mat_list(PyObject *self, PyObject *args){
	unsigned uns;
	PMatArrayObj *pmat_obj = accessor_args_to_PMat(self, args);
	if (pmat_obj == 0L)
		return 0L;
	uns = pmat_obj->n_states;
	return double3DMatToList((const double ***)pmat_obj->p_mat, pmat_obj->n_mat, uns, uns);
}
static PyObject* cpmat_array_get_mat(PyObject *self, PyObject *args){
	unsigned uns;
	int i;
	PyObject *pmat_array_py_obj;
	PMatArrayObj *pmat_array_obj;
	if (!PyArg_ParseTuple(args, "O!i", &pmat_array_type, &pmat_array_py_obj, &i))
		return 0L;
	pmat_array_obj = (PMatArrayObj *)(pmat_array_py_obj);
	i = pyIndToCInd(i, pmat_array_obj->n_mat);
	if (i < 0) {
		return 0L;
	}
	uns = pmat_array_obj->n_states;
	return doubleMatToList((const double **)pmat_array_obj->p_mat[i], uns, uns);
}
/* the factory function
*/
static PyObject* casrvo_ctor(PyObject *self, PyObject *args) {
	int n_categ, mode_enum;
	double val;
	unsigned dim;
	if (!PyArg_ParseTuple(args, "dii", &val, &n_categ, &mode_enum))
		return 0;
	if (n_categ < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of categories parameter must be greater than 0");
		return 0;
	}
	dim = (unsigned) n_categ;
	return (PyObject*) asrv_obj_new(dim, mode_enum, val);
}
static PyObject* casrvo_get_shape(PyObject *self, PyObject *args) {
	PyObject *asrh_py_obj;
	ASRVObj * asrv;
	if (!PyArg_ParseTuple(args, "O!", &asrv_type, &asrh_py_obj))
		return 0L;
	asrv = (ASRVObj *)(asrh_py_obj);
	return PyFloat_FromDouble(asrv->param);
}
static PyObject* casrvo_get_n_cat(PyObject *self, PyObject *args) {
	PyObject *asrh_py_obj;
	ASRVObj * asrv;
	if (!PyArg_ParseTuple(args, "O!", &asrv_type, &asrh_py_obj))
		return 0L;
	asrv = (ASRVObj *)(asrh_py_obj);
	return PyInt_FromLong((long)(asrv->n));
}
static PyObject* casrvo_set_shape(PyObject *self, PyObject *args) {
	double val;
	PyObject *mfo_py_obj;
	if (!PyArg_ParseTuple(args, "O!d", &asrv_type, &mfo_py_obj, &val))
		return 0L;
	if (val <= 0.0) {
		PyErr_SetString(PyExc_ValueError, "The shape parameter must be > 0.0");
		return 0L;
	}
	internal_asrv_set_shape((ASRVObj *)(mfo_py_obj), val);
	return none();
}
static PyObject* casrvo_get_rates(PyObject *self, PyObject *args) {
	PyObject *asrh_py_obj;
	ASRVObj * asrv;
	if (!PyArg_ParseTuple(args, "O!", &asrv_type, &asrh_py_obj))
		return 0L;
	asrv = (ASRVObj *)(asrh_py_obj);
	return doubleArrayToList(asrv->val, asrv->n);
}
static PyObject* calc_pmat(PyObject *self, PyObject *args) {
	int mod_ind;
	double br_len;
	PyObject *pmat_array_py_obj;
	PMatArrayObj *pmat_array_obj;
	PyObject * dsct_model_py_obj;
	DSCTModelObj *dsct_model_obj;
	if (!PyArg_ParseTuple(args, "O!iO!d", &pmat_array_type, &pmat_array_py_obj, &mod_ind, &dsct_model_type, &dsct_model_py_obj, &br_len))
		return 0L;
	pmat_array_obj = (PMatArrayObj *)(pmat_array_py_obj);
	mod_ind = pyIndToCInd(mod_ind, pmat_array_obj->n_mat);
	if (mod_ind < 0) {
		return 0L;
	}
	dsct_model_obj = (DSCTModelObj *)(dsct_model_py_obj);
	pmat_array_obj->model_aliases[mod_ind] = dsct_model_obj;
	pmat_array_obj->brlen_aliases[mod_ind] = br_len;
	if (!do_pmat_calc(pmat_array_obj->p_mat[mod_ind], dsct_model_obj, br_len))
		return 0;
	pmat_array_obj->calc_time_stamp = next_calc_stamp();
	return none();
}

static PyObject* calc_pmat_array(PyObject *self, PyObject *args) {
	unsigned dim, i;
	PyObject * tmp_py_obj;
	PyObject * f_item;
	PyObject * br_len_list_obj = 0L;
	PyObject * model_list_obj = 0L;
	PyObject *pmat_array_py_obj;
	PMatArrayObj *pmat_array_obj;
	if (!PyArg_ParseTuple(args, "O!O!O!", &pmat_array_type, &pmat_array_py_obj, &PyList_Type, &model_list_obj, &PyList_Type, &br_len_list_obj))
		return 0L;
	pmat_array_obj = (PMatArrayObj *)(pmat_array_py_obj);
	dim = pmat_array_obj->n_mat;
	if ((unsigned) PyList_Size(model_list_obj) < dim) {
		PyErr_SetString(PyExc_ValueError, "The list of models must be at least as long as the number of P-matrices");
		return 0L;
	}
	if ((unsigned) PyList_Size(br_len_list_obj) < dim) {
		PyErr_SetString(PyExc_ValueError, "The list of branch lengths must be at least as long as the number of P-matrices");
		return 0L;
	}
	/*Populate the model_aliases and brlen_aliases of the PMatArray*/
	for (i = 0; i < dim; ++i) {
		tmp_py_obj = PyList_GetItem(model_list_obj, i);
		if (tmp_py_obj == 0L) {
			return 0L;
		}
		if (!PyType_IsSubtype(tmp_py_obj->ob_type, &dsct_model_type)){
			PyErr_SetString(PyExc_TypeError, "Expecting a list of DSCTModelObj objects as the second argument");
			return 0L;
		}
		pmat_array_obj->model_aliases[i] = (DSCTModelObj *)tmp_py_obj;
		tmp_py_obj = PyList_GetItem(br_len_list_obj, i);
		if (tmp_py_obj == 0L) {
			return 0L;
		}
		Py_INCREF(tmp_py_obj);
		f_item = PyNumber_Float(tmp_py_obj);
		Py_DECREF(tmp_py_obj);
		if (f_item == 0L) {
			return 0L;
		}
		pmat_array_obj->brlen_aliases[i] = PyFloat_AsDouble(tmp_py_obj);
		Py_DECREF(f_item);
	}
	if (do_pmat_array_calc(pmat_array_obj, 0, dim) == 0L)
		return 0L;
	pmat_array_obj->calc_time_stamp = next_calc_stamp();
	return none();
}
int cla_leaf_consistency_check(CLAObj *c, LeafDataObj *l) {
	assert(c);
	assert(l);
	/* removed this check because we no longer store n_categ at the leaves
	if (c->n_categ < l->n_categ) {
		PyErr_SetString(PyExc_ValueError, "Internal CLA and tip data must have the same n_categ");
		return 0;
	}
	*/
	if (c->n_states < l->n_states) {
		PyErr_SetString(PyExc_ValueError, "Internal CLA and tip data must have the same n_states");
		return 0;
	}
	return 1;
}
int cla_cla_consistency_check(CLAObj *c, CLAObj *o) {
	assert(c);
	assert(o);
	if (c->cso.total_n_sites < o->cso.total_n_sites) {
		PyErr_SetString(PyExc_ValueError, "Internal CLA and tip data must have the same n_sites");
		return 0;
	}
	if (c->cso.n_categ_or_subs < o->cso.n_categ_or_subs) {
		PyErr_SetString(PyExc_ValueError, "Internal CLA and tip data must have the same n_categ");
		return 0;
	}
	if (c->n_states < o->n_states) {
		PyErr_SetString(PyExc_ValueError, "Internal CLA and tip data must have the same n_states");
		return 0;
	}
	return 1;
}

static PyObject* two_tip_cla(PyObject *self, PyObject *args) {
	int ibeg_cat, iend_cat, irescale_thresh, iend_subset, ibeg_subset;
	PyObject *par_cla_py_obj;
	PyObject *f_leaf_data_py_obj;
	PyObject *f_pmat_array_py_obj;
	PyObject *s_leaf_data_py_obj;
	PyObject *s_pmat_array_py_obj;
	CalcContext context;
	CLAObj * par_cla_obj;
	LeafDataObj *f_leaf_data_obj;
	PMatArrayObj *f_pmat_array_obj;
	LeafDataObj *s_leaf_data_obj;
	PMatArrayObj *s_pmat_array_obj;
	int rc, n_subsets;
	if (!PyArg_ParseTuple(args, "O!O!O!O!O!iiiii", &cla_type, &par_cla_py_obj, &leaf_data_type, &f_leaf_data_py_obj, &pmat_array_type, &f_pmat_array_py_obj, &leaf_data_type, &s_leaf_data_py_obj, &pmat_array_type, &s_pmat_array_py_obj, &ibeg_subset, &iend_subset, &ibeg_cat, &iend_cat, &irescale_thresh))
		return 0L;
	par_cla_obj = (CLAObj *)(par_cla_py_obj);
	n_subsets = (par_cla_obj->cso.n_categ_arr == 0L ? 1 : par_cla_obj->cso.n_categ_or_subs);
	f_leaf_data_obj = (LeafDataObj *)(f_leaf_data_py_obj);
	f_pmat_array_obj = (PMatArrayObj *)(f_pmat_array_py_obj);
	s_leaf_data_obj = (LeafDataObj *)(s_leaf_data_py_obj);
	s_pmat_array_obj = (PMatArrayObj *)(s_pmat_array_py_obj);
	if (!cla_leaf_consistency_check(par_cla_obj, f_leaf_data_obj))
		return 0L;
	if (!cla_leaf_consistency_check(par_cla_obj, s_leaf_data_obj))
		return 0L;
	/*@TMP@ should not be using max_n_categ here when the number of subsets is not 1 and there is variability in the # of categories per subset. */
	if (!configure_context(ibeg_subset, iend_subset, n_subsets, ibeg_cat, iend_cat, irescale_thresh, par_cla_obj->cso.max_n_categ, &context))
		return 0L;
	rc = do_two_tip_cla(f_leaf_data_obj, f_pmat_array_obj, s_leaf_data_obj, s_pmat_array_obj, par_cla_obj, context);
#	if defined(PRINTING_LOTS) && PRINTING_LOTS
		print_cla(par_cla_obj);
#	endif	
	if (rc == UNDERFLOW_ERROR_CODE) {
		PyErr_SetString(CLAUnderflowError, "Underflow detected.");
		return 0L;
	}
	return PyInt_FromLong((long)irescale_thresh);
}

static PyObject* one_tip_one_internal_cla(PyObject *self, PyObject *args) {
	int ibeg_cat, iend_cat, irescale_thresh, iend_subset, ibeg_subset;
	PyObject *par_cla_py_obj;
	PyObject *leaf_data_py_obj;
	PyObject *leaf_pmat_array_py_obj;
	PyObject *child_cla_py_obj;
	PyObject *pmat_array_py_obj;
	CalcContext context;
	CLAObj * par_cla_obj;
	LeafDataObj *leaf_data_obj;
	PMatArrayObj *leaf_pmat_array_obj;
	CLAObj *child_cla_obj;
	PMatArrayObj *pmat_array_obj;
	int rc, n_subsets;
	if (!PyArg_ParseTuple(args, "O!O!O!O!O!iiiii", &cla_type, &par_cla_py_obj, &leaf_data_type, &leaf_data_py_obj, &pmat_array_type, &leaf_pmat_array_py_obj, &cla_type, &child_cla_py_obj, &pmat_array_type, &pmat_array_py_obj, &ibeg_subset, &iend_subset, &ibeg_cat, &iend_cat, &irescale_thresh))
		return 0L;
	par_cla_obj = (CLAObj *)(par_cla_py_obj);
	n_subsets = (par_cla_obj->cso.n_categ_arr == 0L ? 1 : par_cla_obj->cso.n_categ_or_subs);
	leaf_data_obj = (LeafDataObj *)(leaf_data_py_obj);
	leaf_pmat_array_obj = (PMatArrayObj *)(leaf_pmat_array_py_obj);
	child_cla_obj = (CLAObj *)(child_cla_py_obj);
	pmat_array_obj = (PMatArrayObj *)(pmat_array_py_obj);
	if (!cla_leaf_consistency_check(par_cla_obj, leaf_data_obj))
		return 0L;
	if (!cla_cla_consistency_check(par_cla_obj, child_cla_obj))
		return 0L;
	/*@TMP@ should not be using max_n_categ here when the number of subsets is not 1 and there is variability in the # of categories per subset. */
	if (!configure_context(ibeg_subset, iend_subset, n_subsets, ibeg_cat, iend_cat, irescale_thresh, par_cla_obj->cso.max_n_categ, &context))
		return 0L;
	rc = do_one_tip_cla(child_cla_obj, pmat_array_obj, leaf_data_obj, leaf_pmat_array_obj, par_cla_obj, context);
	if (rc == UNDERFLOW_ERROR_CODE) {
		PyErr_SetString(CLAUnderflowError, "Underflow detected.");
		return 0L;
	}
	return PyInt_FromLong((long)irescale_thresh);
}
static PyObject* two_internal_cla(PyObject *self, PyObject *args) {
	int ibeg_cat, iend_cat, irescale_thresh, iend_subset, ibeg_subset;
	PyObject *par_cla_py_obj;
	PyObject *f_cla_py_obj;
	PyObject *f_pmat_array_py_obj;
	PyObject *s_cla_py_obj;
	PyObject *s_pmat_array_py_obj;
	CalcContext context;
	CLAObj * par_cla_obj;
	CLAObj * f_cla_obj;
	PMatArrayObj * f_pmat_array_obj;
	CLAObj * s_cla_obj;
	PMatArrayObj * s_pmat_array_obj;
	int rc, n_subsets;
	if (!PyArg_ParseTuple(args, "O!O!O!O!O!iiiii", &cla_type, &par_cla_py_obj, &cla_type, &f_cla_py_obj, &pmat_array_type, &f_pmat_array_py_obj, &cla_type, &s_cla_py_obj, &pmat_array_type, &s_pmat_array_py_obj, &ibeg_subset, &iend_subset,  &ibeg_cat, &iend_cat, &irescale_thresh))
		return 0L;
	par_cla_obj = (CLAObj *)(par_cla_py_obj);
	n_subsets = (par_cla_obj->cso.n_categ_arr == 0L ? 1 : par_cla_obj->cso.n_categ_or_subs);
	f_cla_obj = (CLAObj *)(f_cla_py_obj);
	f_pmat_array_obj = (PMatArrayObj *)(f_pmat_array_py_obj);
	s_cla_obj = (CLAObj *)(s_cla_py_obj);
	s_pmat_array_obj = (PMatArrayObj *)(s_pmat_array_py_obj);
	if (!cla_cla_consistency_check(par_cla_obj, f_cla_obj))
		return 0L;
	if (!cla_cla_consistency_check(par_cla_obj, s_cla_obj))
		return 0L;
	/*@TMP@ should not be using max_n_categ here when the number of subsets is not 1 and there is variability in the # of categories per subset. */
	if (!configure_context(ibeg_subset, iend_subset, n_subsets, ibeg_cat, iend_cat, irescale_thresh, par_cla_obj->cso.max_n_categ,&context))
		return 0L;
	rc = do_internal_cla(f_cla_obj, f_pmat_array_obj, s_cla_obj, s_pmat_array_obj, par_cla_obj, context);
	if (rc == UNDERFLOW_ERROR_CODE) {
		PyErr_SetString(CLAUnderflowError, "Underflow detected.");
		return 0L;
	}
	return PyInt_FromLong((long)irescale_thresh);
}

static PyObject* add_leaf_to_cla(PyObject *self, PyObject *args) {
	int ibeg_cat, iend_cat, irescale_thresh, iend_subset, ibeg_subset;
	PyObject *par_cla_py_obj;
	PyObject *leaf_data_py_obj;
	PyObject *leaf_pmat_array_py_obj;
	CalcContext context;
	CLAObj * par_cla_obj;
	LeafDataObj *leaf_data_obj;
	PMatArrayObj *leaf_pmat_array_obj;
	int rc, n_subsets;
	if (!PyArg_ParseTuple(args, "O!O!O!iiiii", &cla_type, &par_cla_py_obj, &leaf_data_type, &leaf_data_py_obj, &pmat_array_type, &leaf_pmat_array_py_obj, &ibeg_subset, &iend_subset, &ibeg_cat, &iend_cat, &irescale_thresh))
		return 0L;
	par_cla_obj = (CLAObj *)(par_cla_py_obj);
	n_subsets = (par_cla_obj->cso.n_categ_arr == 0L ? 1 : par_cla_obj->cso.n_categ_or_subs);
	leaf_data_obj = (LeafDataObj *)(leaf_data_py_obj);
	leaf_pmat_array_obj = (PMatArrayObj *)(leaf_pmat_array_py_obj);
	if (!cla_leaf_consistency_check(par_cla_obj, leaf_data_obj))
		return 0L;
	/*@TMP@ should not be using max_n_categ here when the number of subsets is not 1 and there is variability in the # of categories per subset. */
	if (!configure_context(ibeg_subset, iend_subset, n_subsets, ibeg_cat, iend_cat, irescale_thresh, par_cla_obj->cso.max_n_categ, &context))
		return 0L;
	rc = do_add_leaf_to_cla(leaf_data_obj, leaf_pmat_array_obj, par_cla_obj, context);
	if (rc == UNDERFLOW_ERROR_CODE) {
		PyErr_SetString(CLAUnderflowError, "Underflow detected.");
		return 0L;
	}
	return PyInt_FromLong((long)irescale_thresh);
}

static PyObject* add_internal_to_cla(PyObject *self, PyObject *args) {
	int ibeg_cat, iend_cat, irescale_thresh, ibeg_subset, iend_subset;
	PyObject *par_cla_py_obj;
	PyObject *child_cla_py_obj;
	PyObject *pmat_array_py_obj;
	CalcContext context;
	CLAObj * par_cla_obj;
	CLAObj *child_cla_obj;
	PMatArrayObj *pmat_array_obj;
	int rc, n_subsets;
	if (!PyArg_ParseTuple(args, "O!O!O!iiiii", &cla_type, &par_cla_py_obj, &cla_type, &child_cla_py_obj, &pmat_array_type, &pmat_array_py_obj, &ibeg_subset, &iend_subset, &ibeg_cat, &iend_cat, &irescale_thresh))
		return 0L;
	par_cla_obj = (CLAObj *)(par_cla_py_obj);
	n_subsets = (par_cla_obj->cso.n_categ_arr == 0L ? 1 : par_cla_obj->cso.n_categ_or_subs);
	child_cla_obj = (CLAObj *)(child_cla_py_obj);
	pmat_array_obj = (PMatArrayObj *)(pmat_array_py_obj);
	if (!cla_cla_consistency_check(par_cla_obj, child_cla_obj))
		return 0L;
	/*@TMP@ should not be using max_n_categ here when the number of subsets is not 1 and there is variability in the # of categories per subset. */
	if (!configure_context(ibeg_subset, iend_subset, n_subsets, ibeg_cat, iend_cat, irescale_thresh, par_cla_obj->cso.max_n_categ, &context))
		return 0L;
	rc = do_add_internal_to_cla(child_cla_obj, pmat_array_obj, par_cla_obj, context);
	if (rc == UNDERFLOW_ERROR_CODE) {
		PyErr_SetString(CLAUnderflowError, "Underflow detected.");
		return 0L;
	}
	return PyInt_FromLong((long)irescale_thresh);
}


static PyObject* get_cla_vals(PyObject *self, PyObject *args) {
	unsigned n_sites, n_categ, n_states, i, state, len;
	PyObject *par_cla_py_obj;
	PyObject *cla_list;
	PyObject *rescaler_list;
	PyObject *el_obj;
	PyObject *to_return;
	const cla_float_t * arr;
#	if SEPARATE_RESCALER_ARRAY
		rescale_history_t * rescalings;
#	endif
	CLAObj * par_cla_obj;
	if (!PyArg_ParseTuple(args, "O!", &cla_type, &par_cla_py_obj))
		return 0L;
	par_cla_obj = (CLAObj *) par_cla_py_obj;
	if (par_cla_obj->cso.n_categ_arr) {
		PyErr_SetString(PyExc_ValueError, "get_cla_vals is not supported when using paritioned models");
		return 0L;
	}

	n_sites = par_cla_obj->cso.total_n_sites;
	n_categ = par_cla_obj->cso.n_categ_or_subs;
	n_states = par_cla_obj->n_states;
	len = n_sites*n_categ;
	arr = (const cla_float_t *) (par_cla_obj->cla);
	cla_list = PyList_New(len*n_states);
	if (cla_list == 0L)
		return 0L;
	rescaler_list = PyList_New(n_sites*n_categ);
	if (rescaler_list == 0L) {
		Py_DECREF(cla_list);
		return 0L;
	}
	for (i = 0; i < len; ++i) {
		for (state = 0; state < n_states; ++state) {
			el_obj = PyFloat_FromDouble(*arr++);
			if (el_obj == 0L) {
				Py_DECREF(cla_list);
				Py_DECREF(rescaler_list);
				return 0L;
			}
			PyList_SetItem(cla_list, n_states*i + state, el_obj);
		}
#		if ! SEPARATE_RESCALER_ARRAY
			el_obj = PyFloat_FromDouble(*arr++);
			if (el_obj == 0L) {
				Py_DECREF(cla_list);
				Py_DECREF(rescaler_list);
				return 0L;
			}
			PyList_SetItem(rescaler_list, i, el_obj);
#		endif
	}
#	if SEPARATE_RESCALER_ARRAY
		rescalings = cla->rescalings;
		for (i = 0; i < len; ++i) {
			if (rescalings) {
				el_obj = PyInt_FromLong((long)(*rescalings++));
			}
			else {
				el_obj = PyInt_FromLong(0L);
			}
			if (el_obj == 0L) {
				Py_DECREF(cla_list);
				Py_DECREF(rescaler_list);
				return 0L;
			}
			PyList_SetItem(rescaler_list, i, el_obj);
		}
#	endif
	to_return = PyTuple_New(2); 
	PyTuple_SetItem(to_return, 0, cla_list); 
	PyTuple_SetItem(to_return, 1, rescaler_list); 
	return to_return;
}
static PyObject* get_full_la_vals(PyObject *self, PyObject *args) {
	unsigned n_sites, n_categ, n_states, i, state, len;
	double categ_freq;
	const double ** freq_ptr;
	PyObject *full_la_py_obj;
	PyObject *cla_list;
	PyObject *rescaler_list;
	PyObject *freq_list;
	PyObject *el_obj;
	PyObject *to_return;
	PyObject *lnL_list;
	const cla_float_t * arr;
#	if SEPARATE_RESCALER_ARRAY
		rescale_history_t * rescalings;
#	endif
	CLAObj * cla_obj;
	FullLAObj * full_la_obj;
	if (!PyArg_ParseTuple(args, "O!", &full_la_type, &full_la_py_obj))
		return 0L;
	full_la_obj = (FullLAObj *) full_la_py_obj;
	cla_obj = full_la_obj->full_cla;
	freq_ptr = (const double **)full_la_obj->state_categ_freqs;
	if (full_la_obj->cso.n_categ_arr) {
		PyErr_SetString(PyExc_ValueError, "get_full_la_vals is not supported when using paritioned models");
		return 0L;
	}
	n_sites = full_la_obj->cso.total_n_sites;
	n_categ = full_la_obj->cso.n_categ_or_subs;
	n_states = full_la_obj->n_states;
	len = n_sites*n_categ;
	arr = (const cla_float_t *) (cla_obj->cla);
	cla_list = PyList_New(len*n_states);
	if (cla_list == 0L)
		return 0L;
	rescaler_list = PyList_New(n_sites*n_categ);
	if (rescaler_list == 0L) {
		Py_DECREF(cla_list);
		return 0L;
	}
	for (i = 0; i < len; ++i) {
		for (state = 0; state < n_states; ++state) {
			el_obj = PyFloat_FromDouble(*arr++);
			if (el_obj == 0L) {
				Py_DECREF(cla_list);
				Py_DECREF(rescaler_list);
				return 0L;
			}
			PyList_SetItem(cla_list, n_states*i + state, el_obj);
		}
#		if ! SEPARATE_RESCALER_ARRAY
			el_obj = PyFloat_FromDouble(*arr++);
			if (el_obj == 0L) {
				Py_DECREF(cla_list);
				Py_DECREF(rescaler_list);
				return 0L;
			}
			PyList_SetItem(rescaler_list, i, el_obj);
#		endif
	}
#	if SEPARATE_RESCALER_ARRAY
		rescalings = cla->rescalings;
		for (i = 0; i < len; ++i) {
			if (rescalings) {
				el_obj = PyInt_FromLong((long)(*rescalings++));
			}
			else {
				el_obj = PyInt_FromLong(0L);
			}
			if (el_obj == 0L) {
				Py_DECREF(cla_list);
				Py_DECREF(rescaler_list);
				return 0L;
			}
			PyList_SetItem(rescaler_list, i, el_obj);
		}
#	endif
	freq_list = PyList_New(n_categ*n_states);
	if (freq_list == 0L) {
		Py_DECREF(cla_list);
		Py_DECREF(rescaler_list);
		return 0L;
	}
	for (i = 0; i < n_categ; ++i) {
		categ_freq = freq_ptr[i][n_states];
		for (state = 0; state < n_states; ++state) {
			el_obj = PyFloat_FromDouble(categ_freq*freq_ptr[i][state]);
			if (el_obj == 0L) {
				Py_DECREF(freq_list);
				Py_DECREF(cla_list);
				Py_DECREF(rescaler_list);
				return 0L;
			}
			PyList_SetItem(freq_list, n_states*i + state, el_obj);
		}	
	}
	arr = (const cla_float_t *) (full_la_obj->pat_lnL_wts);
	lnL_list = PyList_New(n_sites);
	if (lnL_list == 0L) {
		Py_DECREF(freq_list);
		Py_DECREF(cla_list);
		Py_DECREF(rescaler_list);
		return 0L;
	}
	for (i = 0; i < n_sites; ++i) {
		el_obj = PyFloat_FromDouble(arr[2*i]);
		if (el_obj == 0L) {
			Py_DECREF(freq_list);
			Py_DECREF(cla_list);
			Py_DECREF(lnL_list);
			Py_DECREF(rescaler_list);
			return 0L;
		}
		PyList_SetItem(lnL_list, i, el_obj);
	}
	to_return = PyTuple_New(4); 
	PyTuple_SetItem(to_return, 0, lnL_list); 
	PyTuple_SetItem(to_return, 1, freq_list); 
	PyTuple_SetItem(to_return, 2, cla_list); 
	PyTuple_SetItem(to_return, 3, rescaler_list); 
	return to_return;
}
static PyObject* ln_L_internal_edge(PyObject *self, PyObject *args) {
	int ibeg_cat, iend_cat, irescale_thresh, iend_subset, ibeg_subset, n_subsets;
	PyObject *full_la_py_obj;
	PyObject *par_cla_py_obj;
	PyObject *child_cla_py_obj;
	PyObject *pmat_array_py_obj;
	CalcContext context;
	FullLAObj * full_la_obj;
	CLAObj * par_cla_obj;
	CLAObj *child_cla_obj;
	PMatArrayObj *pmat_array_obj;
	double lnL;
	if (!PyArg_ParseTuple(args, "O!O!O!O!iiiii", &full_la_type, &full_la_py_obj,  &cla_type, &child_cla_py_obj, &pmat_array_type, &pmat_array_py_obj, &cla_type, &par_cla_py_obj, &ibeg_subset, &iend_subset, &ibeg_cat, &iend_cat, &irescale_thresh))
		return 0L;
	full_la_obj = (FullLAObj *)(full_la_py_obj);
	par_cla_obj = (CLAObj *)(par_cla_py_obj);
	n_subsets = (par_cla_obj->cso.n_categ_arr == 0L ? 1 : par_cla_obj->cso.n_categ_or_subs);
	child_cla_obj = (CLAObj *)(child_cla_py_obj);
	pmat_array_obj = (PMatArrayObj *)(pmat_array_py_obj);
	if (!cla_cla_consistency_check(par_cla_obj, child_cla_obj))
		return 0L;
	if (!cla_cla_consistency_check(par_cla_obj, full_la_obj->full_cla))
		return 0L;
	/*@TMP@ should not be using max_n_categ here when the number of subsets is not 1 and there is variability in the # of categories per subset. */
	if (!configure_context(ibeg_subset, iend_subset, n_subsets, ibeg_cat, iend_cat, irescale_thresh, par_cla_obj->cso.max_n_categ, &context))
		return 0L;
	lnL = internal_edge_ln_likelihood(child_cla_obj, pmat_array_obj, par_cla_obj, full_la_obj, context);
	return PyFloat_FromDouble(lnL);
}

static PyObject* ln_L_terminal_edge(PyObject *self, PyObject *args) {
	int ibeg_cat, iend_cat, irescale_thresh, iend_subset, ibeg_subset, n_subsets;
	PyObject *full_la_py_obj;
	PyObject *par_cla_py_obj;
	PyObject *tip_data_py_obj;
	PyObject *pmat_array_py_obj;
	CalcContext context;
	FullLAObj * full_la_obj;
	CLAObj * par_cla_obj;
	LeafDataObj *leaf_data_obj;
	PMatArrayObj *pmat_array_obj;
	double lnL;
	if (!PyArg_ParseTuple(args, "O!O!O!O!iiiii", &full_la_type, &full_la_py_obj,  &leaf_data_type, &tip_data_py_obj, &pmat_array_type, &pmat_array_py_obj, &cla_type, &par_cla_py_obj, &ibeg_subset, &iend_subset, &ibeg_cat, &iend_cat, &irescale_thresh))
		return 0L;
	full_la_obj = (FullLAObj *)(full_la_py_obj);
	par_cla_obj = (CLAObj *)(par_cla_py_obj);
	n_subsets = (par_cla_obj->cso.n_categ_arr == 0L ? 1 : par_cla_obj->cso.n_categ_or_subs);
	leaf_data_obj = (LeafDataObj *)(tip_data_py_obj);
	pmat_array_obj = (PMatArrayObj *)(pmat_array_py_obj);
	if (!cla_leaf_consistency_check(par_cla_obj, leaf_data_obj))
		return 0L;
	if (!cla_cla_consistency_check(par_cla_obj, full_la_obj->full_cla))
		return 0L;
	/*@TMP@ should not be using max_n_categ here when the number of subsets is not 1 and there is variability in the # of categories per subset. */
	if (!configure_context(ibeg_subset, iend_subset, n_subsets, ibeg_cat, iend_cat, irescale_thresh, par_cla_obj->cso.max_n_categ, &context))
		return 0L;
	lnL = terminal_edge_ln_likelihood(leaf_data_obj, pmat_array_obj, par_cla_obj, full_la_obj, context);
	return PyFloat_FromDouble(lnL);
}



static PyMethodDef dsct_model_module_functions[] = {
	{"get_full_la_vals", get_full_la_vals, METH_VARARGS,
		"Takes a FullLAObj returns a list of pattern lnL, a list of frequencies for categ x state, the \"root\" conditional likelihoods, and a list of rescalers (either ints or floats).\nCLA's are arranged site x categ x state"},
	{"get_cla_vals", get_cla_vals, METH_VARARGS,
		"Takes a cla object and returns a list of conditional likelihoods, and a list of rescalers (either ints or floats).\nCLA's are arranged site x categ x state"},
	{"calc_ln_L_across_term", ln_L_terminal_edge,	 METH_VARARGS,
		"Calculates ln  likelihood across a terminal branch: args (full LA obj, tip data, tip pmat, par cla, beg cat, end cat, rescale_threshold); Returns ln Likelihood."},
	{"calc_ln_L_across_internal", ln_L_internal_edge,	 METH_VARARGS,
		"Calculates ln  likelihood across an internal branch: args (full LA obj, child cla, pmat, par cla, beg cat, end cat, rescale_threshold); Returns ln Likelihood."},
	{"calc_anc_from_two_tips", two_tip_cla,	 METH_VARARGS,
		"Calculates conditional likelihood array from two tips: args (par cla obj, first tip data, first pmat, second tip data, second pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"calc_anc_from_one_tip_one_intern", one_tip_one_internal_cla,	 METH_VARARGS,
		"Calculates conditional likelihood array from 1 tip and one internal: args (par cla obj,  tip data, tip pmat, child cla, internal pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"calc_anc_from_two_internals", two_internal_cla,	 METH_VARARGS,
		"Calculates conditional likelihood array from two tips: args (par cla obj, first child cla, first pmat, second child cla, second pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"add_leaf_to_cla", add_leaf_to_cla,	 METH_VARARGS,
		"Conditions cla on another tip: args (par cla obj,  tip data,  pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"add_internal_to_cla", add_internal_to_cla,	 METH_VARARGS,
		"Conditions cla on another internal: args (par cla obj,  child cla,  pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"calc_pmat", calc_pmat,	 METH_VARARGS,
		"Calculates a P-Matrix from a Q-matrix and branch length.\nTakes a pmat_array_type, dsct_model_type, and branch length"},
	{"calc_pmat_array", calc_pmat_array,	 METH_VARARGS,
		"Calculates a P-Matrix from a Q-matrix and branch length.\nTakes a pmat_array_type, list of dsct_model_type, and list of branch length"},
	{"csslookup_ctor", csslookup_ctor,	 METH_VARARGS,
		"initializer -- takes a list of lists.  Returns a StateSetLookupStruct."},
	{"cleaf_data_ctor", cleaf_data_ctor, METH_VARARGS,
		"initializer -- takes a list of ints and StateSetLookupStruct"},
	{"cla_ctor", cla_ctor, METH_VARARGS,
		"initializer -- takes n_sites, n_states, n_categ"},
	{"full_la_set_freqs", full_la_set_freqs, METH_VARARGS,
		"Takes a full_la object, a list of category frequencies and a list of lists of state frequencies"},
	{"full_la_ctor", full_la_ctor, METH_VARARGS,
		"initializer -- takes n_sites, n_states, n_categ"},
	{"cdsctm_ctor", cdsctm_ctor,	 METH_VARARGS,
		"initializer -- takes the number of states (integer)."},
	{"cdsctm_get_dimen", cdsctm_get_dimen, METH_VARARGS,
		"returns the number of states"},
	{"cdsctm_get_q_mat", cdsctm_get_q_mat, METH_VARARGS,
		"Returns the Q-matrix as a list of lists of floats"},
	{"cdsctm_get_q_mat_row", cdsctm_get_q_mat_row, METH_VARARGS,
		"Returns a row the Q-matrix as a list of floats (takes row index)."},
	{"cdsctm_get_q_mat_cell", cdsctm_get_q_mat_cell, METH_VARARGS,
		"Returns a cell of the Q-matrix (takes row and col. indices)"},
	{"cdsctm_set_q_mat", cdsctm_set_q_mat, METH_VARARGS,
		"Replaces the Q-matrix (takes a list of list of floats)."},
	{"cdsctm_set_q_mat_row", cdsctm_set_q_mat_row, METH_VARARGS,
		"Replaces a row of the Q-matrix (takes a row index and list)"},
	{"cdsctm_set_q_mat_cell", cdsctm_set_q_mat_cell, METH_VARARGS,
		"Replaces a row of the Q-matrix (takes a row and col. index and value)"},
	{"cpmat_array_ctor", cpmat_array_ctor,	 METH_VARARGS,
		"Initializer -- takes the # of matrices and # of states (positive integers)."},
	{"cpmat_array_get_n_matrices", cpmat_array_get_n_matrices, METH_VARARGS,
		"Returns the "},
	{"cpmat_array_get_n_states", cpmat_array_get_n_states, METH_VARARGS,
		"Returns the number of states"},
	{"cpmat_array_get_mat_list", cpmat_array_get_mat_list, METH_VARARGS,
		"Returns a list P-matrix as a list of  list of lists of floats"},
	{"cpmat_array_get_mat", cpmat_array_get_mat, METH_VARARGS,
		"Takes an index and returns a P-matrix as a list of lists of floats"},
	{"casrvo_ctor", casrvo_ctor, METH_VARARGS,
		"intializer -- takes shape parameter, ncat, dcst_model.RateHetType facet."},
	{"casrvo_get_n_cat", casrvo_get_n_cat, METH_VARARGS,
		"Returns the number of categories."},
	{"casrvo_get_rates", casrvo_get_rates, METH_VARARGS,
		"Returns a list of the rate multipliers."},
	{"casrvo_get_shape", casrvo_get_shape, METH_VARARGS,
		"Returns the shape parameter."},
	{"casrvo_set_shape", casrvo_set_shape, METH_VARARGS,
		"Sets the shape parameter (should be a positive double)."},
/*	{"set_car", set_car, METH_VARARGS},
	{"set_cdr", set_cdr, METH_VARARGS},
*/	{0, 0}
};
/* module entry-point (module-initialization) function */
void
initdsct_model(void)
{
	PyObject * d, * m;
		/* Create the module, with its functions */
	m = Py_InitModule("dsct_model", dsct_model_module_functions);
		/* Finish initializing the type-objects */
	state_set_lookup_type.ob_type = &PyType_Type;
	leaf_data_type.ob_type = &PyType_Type;
	cla_type.ob_type = &PyType_Type;
	full_la_type.ob_type = &PyType_Type;
	dsct_model_type.ob_type = &PyType_Type;
	pmat_array_type.ob_type = &PyType_Type;
	asrv_type.ob_type = &PyType_Type;

	d = PyModule_GetDict(m); 

	CLAUnderflowError = PyErr_NewException("claUnderflow.error", NULL, NULL); 
	PyDict_SetItemString(d, "error", CLAUnderflowError); 

}


/*
################################################################################
# cPhyProb is a package implementing some probability calculations used in
#   calculating likelihoods on phylogenies.
#
# Copyright (C) 2005-2007  Mark Holder mtholder@gmail.com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU  General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
# You should have received a copy of the GNU  General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
################################################################################
*/
