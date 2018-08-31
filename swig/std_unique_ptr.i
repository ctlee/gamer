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

// https://stackoverflow.com/questions/27693812/how-to-handle-unique-ptrs-with-swig
namespace std {
  %feature("novaluewrapper") unique_ptr;
  template <typename Type>
  struct unique_ptr {
    using pointer = Type*;

    explicit unique_ptr( pointer Ptr );
    unique_ptr (unique_ptr&& Right);
    template<class Type2, Class Del2> unique_ptr( unique_ptr<Type2, Del2>&& Right );
    unique_ptr( const unique_ptr& Right) = delete;

    pointer operator-> () const;
    pointer release ();
    void reset (pointer __p=pointer());
    void swap (unique_ptr &__u);
    pointer get () const;
    operator bool () const;

    ~unique_ptr();
  };
}

%define wrap_unique_ptr(Name, Type)
  %template(Name) std::unique_ptr<Type>;
  %newobject std::unique_ptr<Type>::release;  // necessary to prevent memory leak

  %typemap(out) std::unique_ptr<Type> %{
    $result = SWIG_NewPointerObj(new $1_ltype(std::move($1)), $&1_descriptor, SWIG_POINTER_OWN);
  %}
%enddef