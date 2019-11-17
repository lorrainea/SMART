/**
    SMART: Supermaximal Approximate Repeats Tool.
    Copyright (C) 2019 Ayad, L.A.K., Charalampopoulos, P. and Pissis, S.P. 

    ALFRED: A practical method for alignment-free distance computation.
    Copyright (C) 2015 Thankachan, S.V., Chockalingam, S.P., Liu, Y., 
    Apostolico, A. and Aluru, S

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <sstream>
#include "defs.hpp"

using namespace std;

// UTILITY FUNCTIONS -------------------------------------------------
// trim taken from stack overflow
// http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
            s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}


static inline bool ends_with(std::string const &value,
                             std::string const &ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

static inline bool is_valid_dna(char c, std::string& alphabet){
    for(auto i = 0u; i < DNA_ALPHABET_SIZE; i++)
        if(toupper(c) == DNA_ALPHABET[i])
	{
	    if( alphabet == "" )
	    	alphabet = "DNA";
            return true;
	}
    return false;
}

static inline bool is_valid_prot(char c, std::string& alphabet){
    for(auto i = 0u; i < PROT_ALPHABET_SIZE; i++)
        if(toupper(c) == PROT_ALPHABET[i])
	{
            if( alphabet == "DNA" || alphabet == "")
	    	alphabet = "PROT";
            return true;
	}
    return false;
}

static inline bool is_valid_apha(char c, std::string& alphabet){
    return is_valid_dna(toupper(c), alphabet) || is_valid_prot(toupper(c), alphabet);
}

#endif /* UTIL_H */
