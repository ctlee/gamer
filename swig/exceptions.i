/* -*- C -*- */
// Copyright (C) 2012 Johan Hake
//
// This file is part of GAMER.
//
// GAMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GAMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GAMER. If not, see <http://www.gnu.org/licenses/>.

// ===========================================================================
// SWIG directives for exception handling in PyGAMER
// ===========================================================================

// ---------------------------------------------------------------------------
// Function that handles exceptions. Reduces code bloat.
// ---------------------------------------------------------------------------
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
    PyErr_SetString(PyExc_StandardError, const_cast<char*>(e.what()));
  }

  // all runtime_error subclasses
  catch (std::runtime_error &e) 
  {
    PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
  }

  // all the rest
  catch (std::exception &e) 
  {
    PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
  }

}
%}

// ---------------------------------------------------------------------------
// Define the code that each call to GAMER should be wrapped in
// ---------------------------------------------------------------------------
%exception {
  try 
  {
    $action
  }
  catch (...)
  {
    
    // No need to call PyErr_SetString if the error originates from Python
    if (!PyErr_Occurred())
      handle_gamer_exceptions();
    SWIG_fail;
  }
}


