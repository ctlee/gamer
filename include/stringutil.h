#pragma once

#include <algorithm>
#include <cctype>
#include <functional>
#include <string>

namespace stringutil
{
    inline std::vector<std::string> split(const std::string& cstr, std::vector<char> delim = {' ', '\t'})
    {
        std::string s = cstr;
        // convert all delims into delim[0]
        for (auto i=1; i < delim.size(); ++i){
            std::replace(s.begin(), s.end(), delim[i], delim[0]);
        }
        std::vector<std::string> result;
        auto begin = s.begin();
        do{
            auto end = begin;
            while(*end != delim[0] && end != s.end())
                end++;
            if(end != begin) 
                result.push_back(std::string(begin,end));
            begin = end;
        } while (begin++ != s.end());  
        return result;
    }

    std::function<int(int)> isntspace = [](int c) -> int{
        return !std::isspace(c);
    };

    inline void ltrim(std::string& s)
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), isntspace));
    }

    inline void rtrim(std::string& s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(), isntspace).base(), s.end());
    }

    inline void trim(std::string& s)
    {
        ltrim(s);
        rtrim(s);
    }
}