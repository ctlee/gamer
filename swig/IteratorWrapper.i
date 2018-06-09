// ***************************************************************************
// This file is part of the GAMer software.
// Copyright (C) 2016-2017
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


%inline %{
class StopIterator {};

template <typename IT, typename T>
class IteratorWrapper {
public:
  // Constructor
  IteratorWrapper(IT _begin, IT _end): curr(_begin), end(_end) {}

  // Required iterator container function
  IteratorWrapper* __iter__(){
    return this;
  }

  // Iterator function for Python 2.x
  T& next(){
    return __next__();
  }

  // Iterator function for Python 3.x
  T& __next__(){
    if(curr != end){
      T& obj = *curr;
      ++curr;
      return obj;
    }
    else{
      throw StopIterator();
    }
  }

private:
  IT curr;
  IT end;
};
%}

// Python 2.x exception
%exception IteratorWrapper::next {
  try
  {
    $action // calls %extend function next() below
  }
  catch (StopIterator)
  {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}

// Python 3.x exception
%exception IteratorWrapper::__next__ {
  try
  {
    $action // calls %extend function next() below
  }
  catch (StopIterator)
  {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}


%inline %{
template <typename IT, typename T>
class IDIteratorWrapper {
public:
  // Constructor
  IDIteratorWrapper(IT _begin, IT _end): curr(_begin), end(_end) {}

  // ~IDIteratorWrapper(){
  //   while(!allocated.empty()){
  //     T* obj = allocated.back();
  //     allocated.pop_back();
  //     delete obj;
  //   }
  // }

  // Required iterator container function
  IDIteratorWrapper* __iter__(){
    return this;
  }

  // Iterator function for Python 2.x
  T next(){
    return __next__();
  }

  // Iterator function for Python 3.x
  T __next__(){
    if(curr != end){
      T obj = std::move(*curr);
      // allocated.push_back(obj);
      ++curr;
      return std::move(obj);
    }
    else{
      throw StopIterator();
    }
  }

private:
  // std::vector<T*> allocated;
  IT curr;
  IT end;
};
%}

// Python 2.x exception
%exception IDIteratorWrapper::next {
  try
  {
    $action // calls %extend function next() below
  }
  catch (StopIterator)
  {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}

// Python 3.x exception
%exception IDIteratorWrapper::__next__ {
  try
  {
    $action // calls %extend function next() below
  }
  catch (StopIterator)
  {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}