// ***************************************************************************
// This file is part of the GAMer software.
// Copyright (C) 2016-2018
// by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
//    and Michael Holst

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

// ***************************************************************************

%{
SWIGINTERN void handle_gamer_exceptions()
{
  // Re-throw any exception
  try
  {
    throw;
  }
  // all logic_error subclasses
  catch (std::logic_error &e)
  {
    PyErr_SetString(PyExc_SyntaxError, const_cast<char*>(e.what()));
  }
  // all runtime_error subclasses
  catch (std::runtime_error &e)
  {
    PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
  }
  catch (std::out_of_range &e)
  {
    PyErr_SetString(PyExc_IndexError, const_cast<char*>(e.what()));
  }
  // all the rest
  catch (std::exception &e)
  {
    PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
  }
}
%}

// Code that each call to GAMER should be wrapped in
%exception {
  try{
    $action
  }
  catch (...){
    // No need to call PyErr_SetString if the error originates from Python
    if (!PyErr_Occurred()){
      handle_gamer_exceptions();
    }
    SWIG_fail;
  }
}
