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


void print_lcpk(const ReadsDB& rdb, const ivec_t lcpKXY[2][2], const unsigned& k, const unsigned& l, std::ostream& lfs, std::string& alphabet)
{

	const std::string& sx = rdb.getReadById(0);

	int a = 0;
	int b = 0;
	int length = lcpKXY[0][0].size();
	double prob = 0;

	int matrix = 0;

	if( alphabet == "DNA")
	{ 
		matrix = 0;
		prob = 0.25;
	}
	else 
	{
		prob = 0.05;
    		matrix = 1;
	}


	for(int i = 0; i<length; i++ )
	{
		int raw = 0;
		double evalue = 0;
		int err  = k;

		a = max( b-1, lcpKXY[0][1][i] );
		
		if( a == lcpKXY[0][1][i] && b-1 != lcpKXY[0][1][i] )
		{	
	
			if( (unsigned) lcpKXY[0][1][i] >= l )
			{
				int pos1 = i;
				int pos2 = lcpKXY[0][0][i];
				int len = lcpKXY[0][1][i];
	
		
				if(  ( length - len ) - max( pos1, pos2 ) < (signed) k )
					err = ( length - len ) - max( pos1, pos2 );

				if( err == 0 )
					evalue = 0.5 * ( length - len + 1 ) * ( length - len ) * pow( prob, len ) * ( 1 - prob ) + ( length - len ) * pow( prob, l + 1 ) ;
				else evalue = 0.5 * length * ( length - 1 ) * nchoosek( len, err ) * pow( prob, len - err ) * pow( 1-prob, err + 2 );

				for(int i = 0; i<len; i++)
					raw = raw + ( matrix ? pro_delta( sx[pos1 + i], sx[pos2 + i ] ) : nuc_delta( sx[pos1 + i ], sx[pos2 + i ] ) ) ;   

				double identity = (  ( len - k )*1.0 / len*1.0 ) * 100;

				lfs << len<<"\t"<<pos1<<"\t"<<len<<"\t"<<pos2<<"\t"<<err<<"\t"<<fixed<< std::setprecision(3)<<identity<<"\t"<<raw;
				lfs<<"\t"<<scientific<< evalue<<endl;
		  
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
    print_lcpk( rdb, lxy.getkLCP(), cfg.kv, cfg.l, cfg.ofs, cfg.a);
}

void compute_klcp(ReadsDB& rdb, AppConfig& cfg)
{
    klcp_pair_factory<ExactLCPk>( rdb, cfg);
}
