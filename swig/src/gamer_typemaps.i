/* -*- C -*- */
// Copyright (C) 2010 Johan Hake
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-08-06
// Last changed: 2012-03-12
//

//=============================================================================
// General typemaps for GAMer
//=============================================================================

//-----------------------------------------------------------------------------
// The typecheck (unsigned int)
//-----------------------------------------------------------------------------
%typecheck(SWIG_TYPECHECK_INTEGER) unsigned int
{
  $1 = PyInt_Check($input) ? 1 : 0;
}

//-----------------------------------------------------------------------------
// The typemap (unsigned int)
//-----------------------------------------------------------------------------
%typemap(in) unsigned int
{
  if (PyInt_Check($input))
  {
    long tmp = static_cast<long>(PyInt_AsLong($input));
    if (tmp>=0)
      $1 = static_cast<unsigned int>(tmp);
    else
      SWIG_exception(SWIG_TypeError, "expected positive 'int' for argument $argnum");
  }
  else
    SWIG_exception(SWIG_TypeError, "expected positive 'int' for argument $argnum");
}

//%typemap(in, numinputs=0) SurfaceMesh*& (SurfaceMesh* tmp_mesh)
//{
//  $1 = &tmp_mesh;
//}
//
//%typemap(argout) SurfaceMesh*&
//{
//  $result = SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $descriptor(SurfaceMesh*), 1 );
//}
//
//%typemap(in, numinputs=0) SurfaceMesh** (SurfaceMesh* tmp_mesh)
//{
//  $1 = &tmp_mesh;
//}
//
//%typemap(argout) SurfaceMesh**
//{
//  $result = SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $descriptor(SurfaceMesh*), 1 );
//}

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (double* pointlist, int numberofcoordinates)
{
  $1 = PyArray_Check($input) ? 1 : 0;
}

%typemap(in) (double* pointlist, int numberofcoordinates)
{
  if (!PyArray_Check($input))
    SWIG_exception(SWIG_TypeError, "contiguous numpy array of 'd' expected. Make sure that the numpy array is contiguous and uses dtype='d'.");
  PyArrayObject *xa = reinterpret_cast<PyArrayObject*>($input);
  if (!(PyArray_ISCONTIGUOUS(xa) && (PyArray_TYPE(xa) == NPY_DOUBLE)))
    SWIG_exception(SWIG_TypeError, "contiguous numpy array of 'd' expected. Make sure that the numpy array is contiguous and uses dtype='d'.");
  $1 = static_cast<double*>(PyArray_DATA(xa));
  $2 = static_cast<int>(PyArray_DIM(xa,0));
}

%typecheck(SWIG_TYPECHECK_POINTER) (SurfaceMesh** surfmeshes, int num_meshes)
{
  $1 = PyList_Check($input) ? 1 : 0;
}

%typemap (in) (SurfaceMesh** surfmeshes, int num_meshes)
{
  if (!PyList_Check($input))
    SWIG_exception(SWIG_TypeError, "List of SurfaceMeshes expected.");
  
  $2 = PyList_Size($input);
  int res = 0;
  PyObject* py_item = 0;
  void* itemp = 0;
  int newmem = 0;
  SurfaceMesh* sm_ptr = 0;

  // Create Surface mesh pointers
  $1 = new SurfaceMesh*[$2];

  for (int i = 0; i < $2; i++)
  {
    py_item = PyList_GetItem($input, i);
    res = SWIG_ConvertPtr(py_item, (void**)&sm_ptr, $descriptor(SurfaceMesh*), 0);
    if (SWIG_IsOK(res))
      $1[i] = sm_ptr;
    else
      SWIG_exception_fail(SWIG_TypeError, "List of SurfaceMeshes expected.");
  }
}

//-----------------------------------------------------------------------------
// Delete temporary data
//-----------------------------------------------------------------------------
%typemap(freearg) (SurfaceMesh**, int)
{
  delete[] $1;
}
