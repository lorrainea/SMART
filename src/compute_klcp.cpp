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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "ExactLCPk.hpp"
#include "EBLOSUM62.hpp"
#include "EDNAFULL.hpp"


void print64_t_lcpk(const ReadsDB& rdb, const ivec64_t lcpKXY[2][2], const unsigned& k, const unsigned& l, std::ostream& lfs, std::string& alphabet, const unsigned& trim)
{

	const std::string& sx = rdb.getReadById(0);

	int64_t a = 0;
	int64_t b = 0;
	int64_t length = lcpKXY[0][0].size();
	double prob = 0;

	int64_t seq_len = sx.length();

	int64_t matrix = 0;

	if( alphabet == "DNA")
	{ 
		matrix = 0;
		prob = 0.25;
	}
	else 
	{
		prob = 0.04;
    		matrix = 1;
	}


	int64_t * errors  = ( int64_t * ) calloc( ( k ) , sizeof(int64_t ) );
	int64_t len = 0;


	for(int64_t i = 0; i<length; i++ )
	{	
		int64_t raw = 0;
		double evalue = 0;
		int64_t err  = k;

		a = max( b-1, lcpKXY[0][1][i] );

		len = lcpKXY[0][1][i];
		
		if( a == lcpKXY[0][1][i] && b-1 != len )
		{	
		
			if( (unsigned) lcpKXY[0][1][i] >= l )
			{
				int64_t pos1 = i;
				int64_t pos2 = lcpKXY[0][0][i];
				int64_t len = lcpKXY[0][1][i];
				int64_t e = 0;

				int64_t p = 0;
				int64_t c  =0;
				int64_t count = 0;

				
				while( count < len )	
				{
					if( sx[pos1+c] != sx[pos2+c] )
					{
						errors[e] = c;
						e = e + 1;
					}

					if( c < len / 2 )
						c = len - p - 1;
					else 
					{
						p = p + 1;
						c = p ;
					}
					count = count + 1;

				}

				err = e;

				if( err == 0 )
					evalue = 0.5 * ( seq_len - len + 1 ) * ( seq_len - len ) * pow( prob, len ) * ( 1 - prob ) + ( seq_len - len ) * pow( prob, len + 1 ) ;
				else evalue = 0.5 * seq_len * ( seq_len - 1 ) * nchoosek( len, err ) * pow( prob, len - err ) * pow( 1-prob, err + 2 );


				if ( trim == 1 )
				{
					if( err == 1 )
					{
						int64_t mins = min( len - errors[0]- 1, errors[0] );

						int64_t new_pos1 = 0;
						int64_t new_pos2 = 0;
						int64_t new_len = 0;
						int64_t new_err = err;


					
						if( mins == len - errors[0]- 1 )
						{
							new_pos1 = pos1 ;
							new_pos2 = pos2 ;
							new_len = errors[0] ;
							new_err = new_err - 1;

						}
						else if ( mins == errors[0] ) 		
						{		
							new_pos1 = pos1+errors[0]+1;
							new_pos2 = pos2+errors[0]+1;
			
							new_len = len - errors[0] - 1 ;
							new_err = new_err - 1;

						}

						double score= 0;

				
						if( new_err == 0 )
							score = 0.5 * ( seq_len - new_len + 1 ) * ( seq_len - new_len ) * pow( prob, new_len ) * ( 1 - prob ) + ( seq_len - new_len ) * pow( prob, new_len + 1 ) ;
						else score = 0.5 * seq_len * ( seq_len - 1 ) * nchoosek( new_len, new_err ) * pow( prob, new_len - new_err ) * pow( 1-prob, new_err + 2 );
				
					
						if( evalue !=0 && score < evalue)
						{

							pos1 = new_pos1;
							pos2 = new_pos2;
							len = new_len;
							evalue = score;
							err = new_err;
						}
					}
					else if ( err > 1 )
					{	
						int64_t b = 0;
						int64_t orig_pos1 = pos1;
						int64_t orig_pos2 = pos2;
						int64_t orig_len = len;

						while( err >= 1)
						{
							int64_t mins = min( len - errors[b]- 1, errors[b] );
							int64_t new_pos1 = 0;
							int64_t new_pos2 = 0;
							int64_t new_len = 0;
							int64_t new_err = 0;

							if( mins == len - errors[b]- 1 )
							{
								new_pos1 = pos1; 
								new_pos2 = pos2;
								new_len = orig_len - ( orig_len - errors[b] )  - (orig_len - len ) ;
								new_err = err - 1;

							}
							else if ( mins == errors[b] ) 		
							{		
								new_pos1 = orig_pos1+ errors[b]+1;
								new_pos2 = orig_pos2+ errors[b]+1;
			
								new_len = len - errors[b] - 1 ;
								new_err = err - 1;

							}

							double score= 0;

							if( new_err == 0 )
								score = 0.5 * ( seq_len - new_len + 1 ) * ( seq_len - new_len ) * pow( prob, new_len ) * ( 1 - prob ) + ( seq_len - new_len ) * pow( prob, new_len + 1 ) ;
							else score = 0.5 * seq_len * ( seq_len - 1 ) * nchoosek( new_len, new_err ) * pow( prob, new_len - new_err ) * pow( 1-prob, new_err + 2 );


							if( evalue != 0 && score < evalue)
							{	pos1 = new_pos1;
								pos2 = new_pos2;
								len = new_len;
								evalue = score;
								err = new_err;

							}


							b = b+1;

							if( b >= e )
								break;
						
						}

					
					}
				}

				for(int64_t i = 0; i<len; i++)
					raw = raw + ( matrix ? pro_delta( sx[pos1 + i], sx[pos2 + i ] ) : nuc_delta( sx[pos1 + i ], sx[pos2 + i ] ) ) ;   

				double identity = (  ( len - k )*1.0 / len*1.0 ) * 100;

				if( len >= l )
				{
					lfs << len<<"\t"<<pos1<<"\t"<<"F"<<"\t"<<len<<"\t"<<pos2<<"\t"<<err<<"\t"<<fixed<< std::setprecision(3)<<identity<<"\t"<<raw;
					lfs<<"\t"<<scientific<< evalue<<endl;
				}
		  
			}
		}
		
		b = a;	
	}

}

template<typename LCPk>
void klcp_pair_factory( ReadsDB& rdb, AppConfig& cfg){

    const std::string& sx = rdb.getReadById(0);

    LCPk lxy(sx, cfg);
    lxy.compute();
    print64_t_lcpk( rdb, lxy.getkLCP(), cfg.kv, cfg.l, cfg.ofs, cfg.a, cfg.t );
}

void compute_klcp(ReadsDB& rdb, AppConfig& cfg)
{
    klcp_pair_factory<ExactLCPk>( rdb, cfg);
}
