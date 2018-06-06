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
