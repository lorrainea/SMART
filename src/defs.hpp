/**
    SMART: Supermaximal Approximate Repeats Tool.
    Copyright (C) 2019 Ayad, L.A.K., Charalampopoulos, P. and Pissis, S.P.

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

#ifndef TYPE_DEFS_HPP
#define TYPE_DEFS_HPP

#include <cstdint>
#include <vector>

#define PRO_SCORING_MATRIX_SIZE 23
#define NUC_SCORING_MATRIX_SIZE 15
#define pro_delta(a,b) EBLOSUM62_matrix[ BLOSUM[(int)(a)] ][ BLOSUM[(int)(b)] ]
#define nuc_delta(a,b) EDNAFULL_matrix[ EDNA[(int)(a)] ][ EDNA[(int)(b)] ]

#define PROT                    "ARNDCQEGHILKMFPSTWYVBZX*"   //Proteins alphabet
#define DNA                     "ATGC"            //DNA alphabet

// CONSTANTS ---------------------------------------------------------
const unsigned DNA_ALPHABET_SIZE           = 4;
const char DNA_ALPHABET[DNA_ALPHABET_SIZE+1] = {'A', 'C', 'G', 'T', 'N'};
const unsigned ST_ALPHABET_SIZE            = DNA_ALPHABET_SIZE + 1;
const char ST_ALPHABET[ST_ALPHABET_SIZE]   = {'A', 'C', 'G', 'T', '$'};

const unsigned PROT_ALPHABET_SIZE              = 23;
const char PROT_ALPHABET[PROT_ALPHABET_SIZE]   = {'A', 'C', 'D', 'E', 'F', 'G',
                                                  'H', 'I', 'K', 'L', 'M', 'N',
                                                  'P', 'Q', 'R', 'S', 'T', 'V',
                                                  'W', 'Y', 'B', 'Z', 'X'};

extern unsigned int EDNA[];
extern unsigned int BLOSUM[];

// int vectors
typedef std::vector<int32_t> ivec_t;
typedef std::vector<int64_t> ivec64_t;

#endif /* TYPE_DEFS_HPP */
