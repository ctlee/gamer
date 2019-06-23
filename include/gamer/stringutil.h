/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2018
 * by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
 *    and Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * ***************************************************************************
 */

/**
 * @file stringutil.h
 * @brief Various string utilities
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <functional>
#include <string>
#include <vector>


/// Namespace for all things gamer
namespace gamer
{

/// Namespace for string utilities
namespace stringutil
{
/**
 * @brief      Split a string into a vector of substrings
 *
 * @param[in]  cstr   Input string to split
 * @param[in]  delim  Vector of delimiters to split at
 *
 * @return     Vector of substrings split at delimiters
 */
inline std::vector<std::string> split(const std::string &cstr, std::vector<char> delim = {' ', '\t'})
{
    std::string s = cstr;
    // convert all delims into delim[0]
    for (int i = 1; i < delim.size(); ++i)
    {
        std::replace(s.begin(), s.end(), delim[i], delim[0]);
    }
    std::vector<std::string> result;
    auto                     begin = s.begin();
    do
    {
        auto end = begin;
        while (*end != delim[0] && end != s.end())
            end++;
        if (end != begin)
            result.push_back(std::string(begin, end));
        begin = end;
    }
    while (begin++ != s.end());
    return result;
}

/// @cond detail
/// Namespace for internal string utility functions
namespace stringutil_detail
{
/// Functor for negated isspace
std::function<int(int)> isntspace = [](int c) -> int{
        return !std::isspace(c);
    };
}
/// @endcond

/**
 * @brief      Inplace removal of white space from the left side of a string
 *
 * @param      s     String to trim
 */
inline void ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), stringutil_detail::isntspace));
}


/**
 * @brief      Inplace removal of white space from the right side of a string
 *
 * @param      s     String to trim
 */
inline void rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), stringutil_detail::isntspace).base(), s.end());
}

/**
 * @brief      Inplace trimming of white space on both sides of a string
 *
 * @param      s     String to trim
 */
inline void trim(std::string &s)
{
    ltrim(s);
    rtrim(s);
}
} // end namespace stringutil
} // end namespace gamer
