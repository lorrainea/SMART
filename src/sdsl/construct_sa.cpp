#include "construct_sa.hpp"
#include "divsufsort.h"
#include "divsufsort64.h"

//! Construct the Suffix Array for a text.
/*!
 * \param c Text (c-string) to calculate the suffix array. The lex. order is given by the ascii-codes of the characters.
 * \param len Length of the text. *(c+len)=0 and for i<len *(c+len)!=0
 * \param sa Reference to a RandomAccessContainer which will contain the result of the calculation.
 * \pre sa.size() has to be equal to len.
 */
void construct_sa(const unsigned char* text, ivec_t::size_type len, ivec_t& sa)
{
    sa.resize(len);
    divsufsort(text, (int32_t*)(&sa[0]), len);
}

void construct_sa64(const unsigned char* text, ivec_t::size_type len, ivec64_t& sa)
{
    sa.resize(len);
    divsufsort64(text, (int64_t*)(&sa[0]), len);
}
