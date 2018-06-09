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

%ignore index_map::operator[](std::size_t idx);
%inline %{
struct index_map {
  std::vector<int> map;
  std::size_t maxsize;

  index_map(std::size_t size){
    maxsize = size;
    map = std::vector<int>(size);
  }

  void updatemap(int marker, PyObject* incoming){
   if (PyList_Check(incoming)) {
      for(Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
        PyObject *value = PyList_GetItem(incoming, i);
        map[PyLong_AsSize_t(value)] = marker;
      }
    }
    else {
      throw std::logic_error("Passed PyObject pointer was not a list!");
    }
  }

  int operator[](std::size_t idx){
    return map[idx];
  }

  void printMap(){
    for(auto item : map){
      std::cout << item << ", ";
    }
    std::cout << std::endl;
  }

  int __getitem__(int key){
    if (key < 0 || key > map.size()){
      throw std::out_of_range("Index out of bounds");
    }
    return map[key];
  }
};
%}