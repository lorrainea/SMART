/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file construct_sa.hpp
    \brief construct_sa.hpp contains an interface to access suffix array construction algorithms
    \author Simon Gog
*/

#ifndef CONSTRUCT_SA_H
#define CONSTRUCT_SA_H

#include "defs.hpp"
#include <cassert>

void construct_sa(const unsigned char* text, ivec_t::size_type len, ivec_t& sa);
void construct_sa64(const unsigned char* text, ivec_t::size_type len, ivec64_t& sa);

template<typename int_vector>
void construct_isa(const int_vector& sa, int_vector& isa){
    typedef typename int_vector::size_type size_type;
    isa.resize(sa.size());
    for (size_type i=0; i < isa.size(); ++i) {
        isa[ sa[i] ] = i;
    }

}

//! Construct the LCP array for text over byte- or integer-alphabet.
/*!	The algorithm computes the lcp array
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \pre Text and Suffix array exist in the cache. Keys:
 *
 *  \post LCP array
 *
 *  \par Time complexity
 *         \f$ \Order{n} \f$
 *  \par Space complexity
 *         \f$ n (\log \sigma + \log n) \f$ bits
 *  \par Reference
 *         Toru Kasai, Gunho Lee, Hiroki Arimura, Setsuo Arikawa, Kunsoo Park:
 *         Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications.
 *         CPM 2001: 181-192
 */

template<typename int_vector>
void construct_lcp_kasai(const char* text, const int_vector& sa,
                         const int_vector& isa, int_vector& lcp){
    typedef typename int_vector::size_type size_type;
    lcp.resize(sa.size());
    for(size_type i = 0; i < sa.size();i++)
        lcp[i] = sa[i];

    // use Kasai algorithm to compute the lcp values
    for (size_type i=0,j=0,sa_1=0,l=0, n=isa.size(); i < n; ++i) {
        sa_1 =  isa[i]; // = isa[i]
        if (sa_1) {
            j = lcp[sa_1-1];
            if (l) --l;
            assert(i!=j);
            while (text[i+l]==text[j+l]) { // i+l < n and j+l < n are not necessary, since text[n]=0 and text[i]!=0 (i<n) and i!=j
                ++l;
            }
            lcp[ sa_1-1 ] = l; //overwrite sa array with lcp values
        } else {
            l = 0;
            lcp[ n-1 ] = 0;
        }
    }

    for (size_type i=sa.size(); i>1; --i) {
        lcp[i-1] = lcp[i-2];
    }
    lcp[0] = 0;
}


//! Construct the LCP array for text over byte- or integer-alphabet.
/*!	The algorithm computes the lcp array and stores it to disk.
 *  \pre Text and Suffix array .
 *  \post LCP array exist in the cache.
 *  \par Time complexity
 *         \f$ \Order{n} \f$
 *  \par Space complexity
 *         \f$ n( \log \sigma + \log \n ) \f$ bits
 *  \par Reference
 *         Juha Kärkkäinen, Giovanni Manzini, Simon J. Puglisi:
 *         Permuted Longest-Common-Prefix Array.
 *         CPM 2009: 181-192
 */
template<typename int_vector>
void construct_lcp_PHI(const char* text, const int_vector& sa, int_vector& lcp){
    typedef typename int_vector::size_type size_type;
    size_type n = sa.size();
    //	(1) Calculate PHI (stored in array plcp)
    int_vector plcp(n);
    lcp.resize(n);
    for (size_type i=0, sai_1 = 0; i < n; ++i) {
        size_type sai = sa[i];
        plcp[ sai ] = sai_1;
        sai_1 = sai;
    }

    //  (2) Calculate permuted LCP array (text order), called PLCP
    size_type max_l = 0;
    for (size_type i=0, l=0; i < n-1; ++i) {
        size_type phii = plcp[i];
        while (text[i+l] == text[phii+l]) {
            ++l;
        }
        plcp[i] = l;
        if (l) {
            max_l = std::max(max_l, l);
            --l;
        }
    }

    lcp[0] = 0;
    for (size_type i=1; i < n; ++i) {
        size_type sai = sa[i];
        lcp[i] = plcp[sai];
    }

}

#endif /* CONSTRUCT_SA_H */
