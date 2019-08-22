/**
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

#include "compute_klcp.hpp"
#include "construct_sa.hpp"
#include "ExactLCPk.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// Construct SA, ISA, LCP, RMQ data structures for the string
ExactLCPk::ExactLCPk(const std::string& sx, AppConfig& cfg) : m_aCfg(cfg)
{
	m_strXY = sx +"#$";
	m_strLengths[0] = sx.size(); m_strLengths[1] = 1;
	m_strLenPfx[0] = sx.size(); m_strLenPfx[1] = sx.size() + 2;
	m_shiftPos[0] = 0; m_shiftPos[1] = m_strLengths[0] + 1;

	// construct SA, ISA, LCP and RMQ
	construct_sa((const unsigned char*)m_strXY.c_str(), m_strXY.size(), m_gsa);
	construct_isa(m_gsa, m_gisa);
	construct_lcp_kasai(m_strXY.c_str(), m_gsa, m_gisa, m_glcp);

	//construct_lcp_PHI(pxy.c_str(), gsa, glcp);
	m_rangeMinQuery = rmq_support_sparse_table<ivec_t, true, ivec_t>(&m_glcp);
	m_nPass = m_passSizes = 0.0;
}


// Chop prefix of the suffixes corresponding to an internal node of the suffix tree
void ExactLCPk::chopPrefix0(const InternalNode& uNode, std::vector<L1Suffix>& leaves)
{
	if(uNode.m_leftBound == -1 || uNode.m_rightBound == -1)
		return;

	assert(uNode.m_rightBound > uNode.m_leftBound);
	assert(uNode.m_leftBound >= 2);
	assert(uNode.m_rightBound < (int32_t)m_gsa.size());

	leaves.clear();
	leaves.resize(uNode.m_rightBound - uNode.m_leftBound + 1);

	//   collect tuples for each position (going left and right)
	//      (i, i', 0/1) i' = gisa[gsa[i] + d + 1]
	for(int32_t idx = uNode.m_leftBound, i = 0; idx <= uNode.m_rightBound; idx++, i++)
	{
		int32_t spos = m_gsa[idx];
		int32_t rs = spos < m_strLengths[0] ? 0 : 1;
       		int32_t epos = spos + uNode.m_stringDepth + 1;
        	// crossing boundary
        	int32_t esa = epos <= m_strLenPfx[rs] ? m_gisa[epos] : -1;
        	L1Suffix cm(spos, esa, rs);
        	leaves[i] = cm;
    	}
    	//   sort tuples by i'
    	std::sort(leaves.begin(), leaves.end());
}

// Select all suffixes corresponding to internal node of the suffix tree
void ExactLCPk::selectSuffixes0(const InternalNode& uNode, std::vector<L1Suffix>& leaves)
{
	if(uNode.m_leftBound == -1 || uNode.m_rightBound == -1)
		return;

    	assert(uNode.m_rightBound > uNode.m_leftBound);
    	assert(uNode.m_leftBound >= 2);
    	assert(uNode.m_rightBound < (int32_t)m_gsa.size());

    	leaves.clear();
    	leaves.resize(uNode.m_rightBound - uNode.m_leftBound + 1);

    	//   collect tuples for each position (going left and right)
    	//      (i, i', 0/1) i' = gisa[gsa[i] + d + 1]
    	for(int32_t idx = uNode.m_leftBound, i = 0; idx <= uNode.m_rightBound; idx++, i++)
	{
        	int32_t spos = m_gsa[idx];
        	int32_t rs = spos < m_strLengths[0] ? 0 : 1;
        	
		// crossing boundary
        	int32_t esa = spos <= m_strLenPfx[rs] ? idx : -1;

        	L1Suffix cm(spos, esa, rs);
        	leaves[i] = cm;
    	}
}

// eliminate duplicate internal nodes
void ExactLCPk::eliminateDupes(std::vector<InternalNode>& uNodes)
{
	// sort
	std::sort(uNodes.begin(), uNodes.end());
	// move up
	unsigned j = 0;

	for(unsigned i = 1; i < uNodes.size(); i++)
	{
        	if(uNodes[i] == uNodes[j])
			continue;
		j += 1;
		uNodes[j] = uNodes[i];
	}
	uNodes.resize(j+1);
}

int32_t ExactLCPk::ansvLB0(int32_t curLeaf, int32_t ansvBound)
{
	if(curLeaf + 1 >= (int32_t)m_gsa.size()) // reached the end
		curLeaf = (int32_t)m_gsa.size() - 2;
	if(curLeaf <= 2) // First two suffixes comes from the separators
		return 2;
	return (ansvBound < 2) ? 2 : ansvBound;
}


int32_t ExactLCPk::ansvRB0(int32_t curLeaf, int32_t ansvBound)
{
	if(curLeaf <= 2) // First two suffixes comes from the separators
		curLeaf = 2;
	if(curLeaf + 1 >= (int32_t)m_gsa.size()) // reached the end
		return curLeaf;
	return ansvBound;
}

// Select internal nodes of the suffix tree
void ExactLCPk::selectInternalNodes0(std::vector<InternalNode>& uNodes)
{
	uNodes.resize(m_gsa.size() - 2);

	std::vector<int32_t> lbounds, rbounds;
	ansvBounds0<int32_t, int32_t>(m_gsa, m_glcp, lbounds, rbounds);
	assert(lbounds.size() == m_gsa.size());
	assert(rbounds.size() == m_gsa.size());
	// Get all the internal nodes

	for(int32_t i = 2; i < (int32_t)m_gsa.size() - 1; i++)
	{
        	uNodes[i-2].m_leftBound = ansvLB0(i+1, lbounds[i+1]);
        	uNodes[i-2].m_rightBound = ansvRB0(i+1, rbounds[i+1]);
       
        	uNodes[i-2].m_stringDepth = m_glcp[i + 1];
        	uNodes[i-2].m_delta = 0;
    	}

	// Eliminate duplicates
	eliminateDupes(uNodes);
}

// Update pass
// Pass through the leaves to update longest matching prefixes
void ExactLCPk::updateExactLCPk(InternalNode& uNode,std::vector<L1Suffix>& leaves)
{
	#ifdef DEBUG
	m_aCfg.lfs << std::endl;
	for(auto cm: leaves)
	{
		cm.dwriteln(m_aCfg.lfs);
    	}
    	m_aCfg.lfs << std::endl;
	#endif
    
	if(leaves.size() == 0)
		return;
	m_nPass += 1;
	m_passSizes += leaves.size();

	// left -> right pass
	updatePass<UpperBoundCheck, IncrPointer>(0, 1, (int32_t)leaves.size(), uNode, leaves
	#ifdef DEBUG
        , "L->R"
	#endif 
	);
    	// right -> left pass
	updatePass<LowerBoundCheck, DecrPointer>((int32_t)leaves.size() - 1, (int32_t)leaves.size() - 2, -1, uNode, leaves
	#ifdef DEBUG
        , "R->L"
	#endif 
	);

	#ifdef DEBUG
	m_aCfg.lfs << std::endl;
	#endif
}


//  Assuming that the trieLeaves are sorted by m_errSAPos
//  In order to maintain consistency,
//    - we start with curLeaf and its next as the tree 'x'
//    - we grow this tree 'x' until we reach its left end
//
int32_t ExactLCPk::leftBoundK(const std::vector<L1Suffix>& trieLeaves, int32_t curLeaf)
{
	if(curLeaf + 1 >= (int32_t)trieLeaves.size()) // reached the end
		curLeaf = (int32_t) trieLeaves.size() - 2;
	if(curLeaf <= 0)
		return 0;
	if(trieLeaves[curLeaf].m_errSAPos < 2 || trieLeaves[curLeaf + 1].m_errSAPos < 2)
		return curLeaf;
	int32_t lpos, idxLCP = rangeMinLCP(trieLeaves[curLeaf].m_errSAPos, trieLeaves[curLeaf + 1].m_errSAPos), srcLCP = idxLCP;
	while(idxLCP >= srcLCP)
	{
		lpos = curLeaf;
		curLeaf -= 1;
		if(curLeaf < 0 || trieLeaves[curLeaf].m_errSAPos < 2)
			break;
		idxLCP = rangeMinLCP(trieLeaves[curLeaf].m_errSAPos, trieLeaves[curLeaf + 1].m_errSAPos);
	}
	return lpos;
}

//  Assuming that the trieLeaves are sorted by m_errSAPos
//  In order to maintain consistency,
//    - we start with curLeaf and its next as the tree 'x'
//    - we grow this tree 'x' until we reach its right end
//
int32_t ExactLCPk::rightBoundK(const std::vector<L1Suffix>& trieLeaves, int32_t curLeaf)
{
	if(curLeaf < 0)
		curLeaf = 1;
	if(curLeaf + 1 >= (int32_t)trieLeaves.size()) // reached the end
		return (int32_t)trieLeaves.size() - 1;
	// compute rpos (TODO:: verify)
	curLeaf += 1;
	if(trieLeaves[curLeaf - 1].m_errSAPos < 2 || trieLeaves[curLeaf].m_errSAPos < 2)
		return curLeaf;
	int32_t rpos,
	idxLCP = rangeMinLCP(trieLeaves[curLeaf - 1].m_errSAPos, trieLeaves[curLeaf].m_errSAPos), srcLCP = idxLCP;
	while(idxLCP >= srcLCP)
	{
		rpos = curLeaf;
		curLeaf += 1;
		if(curLeaf >= (int32_t)trieLeaves.size() || trieLeaves[curLeaf].m_errSAPos < 2)
			break;
		idxLCP = rangeMinLCP(trieLeaves[curLeaf - 1].m_errSAPos, trieLeaves[curLeaf].m_errSAPos);
	}
return rpos;
}

int32_t ExactLCPk::ansvLBK(const std::vector<L1Suffix>& trieLeaves,
                           int32_t curLeaf, int32_t minIdx,
                           int32_t ansvBound){
  if(curLeaf + 1 >= (int32_t)trieLeaves.size()) // reached the end
    curLeaf = (int32_t) trieLeaves.size() - 2;
  if(curLeaf <= 0)
    return 0;
  return ansvBound < minIdx ? minIdx : ansvBound;
}

int32_t ExactLCPk::ansvRBK(const std::vector<L1Suffix>& trieLeaves,
                           int32_t curLeaf, int32_t ansvBound){
  if(curLeaf < 0)
    curLeaf = 1;
  if(curLeaf + 1 >= (int32_t)trieLeaves.size()) // reached the end
    return (int32_t)trieLeaves.size() - 1;
  return ansvBound;
}

//
// Select internal nodes of the trie resulting from chopping the
//  suffixes (i.e. leaves) corresponding to the internal node
//  from the previous iteration.
void ExactLCPk::selectInternalNodesK(const InternalNode& prevNode, const std::vector<L1Suffix>& leaves,std::vector<InternalNode>& trieNodes)
{
	// Assume input suffixes are sorted (each suffix is a leaf in the trie)
	trieNodes.resize(leaves.size());
	unsigned j = 0;
	std::vector<int32_t> lbounds, rbounds;
	struct RMQ 
	{
		ExactLCPk& rpx;
		RMQ(ExactLCPk& px):rpx(px){}
		inline int32_t operator()(const L1Suffix& x, const L1Suffix& y)
		{
			return rpx.rangeMinLCP(x, y);
        	}
   	 };

	RMQ rminq(*this);
	ansvBoundsK<int32_t, int32_t, RMQ>(leaves, rminq, lbounds, rbounds);
	int32_t minValidIdx = 0;
	for(int32_t i = 0; i < (int32_t)leaves.size(); i++)
	{
		if(leaves[i].m_errSAPos >= 2) // find the one which doesn't reach the end
		break;
		minValidIdx++;
	}
	// for each leaf
	for(int32_t i = 0; i < (int32_t)leaves.size(); i++)
	{
		if(leaves[i].m_errSAPos < 2) // skip the one which reach the end
			continue;
        	// - make an internal node
        	InternalNode uNode;
        	uNode.m_leftBound = ansvLBK(leaves, i, minValidIdx, lbounds[i]);
        	uNode.m_rightBound = ansvRBK(leaves, i, rbounds[i]);
        
        	if(uNode.m_leftBound == uNode.m_rightBound) // skip internal node with 1 leaf
            		continue;
        	// prevNode.depth + 1 + get range min lcp of left and right ends
        	uNode.m_stringDepth = rangeMinLCP(leaves[uNode.m_leftBound].m_errSAPos, leaves[uNode.m_rightBound].m_errSAPos);
        	uNode.m_delta = prevNode.m_delta + prevNode.m_stringDepth + 1;
        	trieNodes[j] = uNode;
        	j++;
    	}
    	trieNodes.resize(j);
    	// remove duplicates of internal nodes
	eliminateDupes(trieNodes);
}

// Chop prefix of the suffixes corresponding to an internal node
//  of the suffix trie
void ExactLCPk::chopPrefixK(const InternalNode& uNode,const std::vector<L1Suffix>& inSuffixes,std::vector<L1Suffix>& outSuffixes)
{
	if(uNode.m_leftBound == -1 || uNode.m_rightBound == -1)
		return;

	assert(uNode.m_rightBound > uNode.m_leftBound);
	outSuffixes.clear();
	outSuffixes.resize(uNode.m_rightBound - uNode.m_leftBound + 1);
	#ifdef DEBUG
	m_aCfg.lfs << std::endl << "   ";
	uNode.dwriteln(m_aCfg.lfs);
	m_aCfg.lfs << std::endl;
	#endif
	//   collect tuples for each position (going left and right)
	//      (c, c', 0/1), where c' = gisa[gsa[c] + d + 1]
	unsigned j = 0;
	for(int32_t idx = uNode.m_leftBound; idx <= uNode.m_rightBound; idx++)
	{
		int32_t epx = inSuffixes[idx].m_errSAPos;
		// ignore that has crossed boundary in prev. level
		if(epx < 0 || epx >= (int32_t)m_gsa.size())
			continue;
		int32_t spos = m_gsa[epx];
		int32_t rs = inSuffixes[idx].m_srcStr;
		int32_t epos = spos + uNode.m_stringDepth + 1 ;
		// crossing boundary
		int32_t esa = epos <= m_strLenPfx[rs] ? m_gisa[epos] : -1;

		L1Suffix cm(spos, esa, rs);
		outSuffixes[j] = cm;
		#ifdef DEBUG
		inSuffixes[idx].dwrite(m_aCfg.lfs);
		cm.dwrite(m_aCfg.lfs);
		m_aCfg.lfs << std::endl;
		#endif
		j++;
	}

	outSuffixes.resize(j);
	//   sort tuples by c'
	std::sort(outSuffixes.begin(), outSuffixes.end());
	#ifdef DEBUG
	m_aCfg.lfs << std::endl;
	#endif
}

// Function to compute LCP with no mismatches
void ExactLCPk::compute0()
{
	assert(m_gsa.size() > 2);
	assert(m_gsa.size() == m_glcp.size()); // an assumption of lca

	// resize the LCP arrays
	for(unsigned i = 0; i < 2; i++)
	{
		m_klcpXY[i][0].resize(m_strLengths[i], 0);
		m_klcpXY[i][1].resize(m_strLengths[i], 0);
	}

	InternalNode uNode;
	uNode.m_leftBound = 2;
	uNode.m_rightBound = m_gsa.size() - 1;
	uNode.m_delta = 0;
	#ifdef DEBUG
    	uNode.dwriteln(m_aCfg.lfs);
	#endif
	//   collect tuples for each position (going left and right)
	//      (i, i', 0/1) , where i' = gisa[gsa[i] + d + 1]
	std::vector<L1Suffix> leaves;
	selectSuffixes0(uNode, leaves);
	uNode.m_stringDepth = -1; // just to let update use the LCP
	// update lcp using sorted tuples using a double pass
	updateExactLCPk(uNode, leaves);
}

// Recursive function to compute LCP_k for a given internal node
// and the corresponding chopped suffixes
void ExactLCPk::computeK(InternalNode& uNode, std::vector<L1Suffix>& uLeaves, int searchLevel)
{
	if(searchLevel == 0)
	{
		// update LCP array
		updateExactLCPk(uNode, uLeaves);
		return;
	}
	std::vector<InternalNode> trieNodes;
	selectInternalNodesK(uNode, uLeaves, trieNodes);
	for(auto nit = trieNodes.begin(); nit != trieNodes.end(); nit++)
	{
		std::vector<L1Suffix> trieLeaves;
		chopPrefixK(*nit, uLeaves, trieLeaves);
		computeK(*nit, trieLeaves, searchLevel - 1);
	}
}

// Entry for LCP_k computation
void ExactLCPk::computeK()
{
	assert(m_gsa.size() > 2);
	assert(m_gsa.size() == m_glcp.size()); // an assumption of lca
	// resize the LCP arrays
	for(unsigned i = 0; i < 2; i++)
	{
		m_klcpXY[i][0].resize(m_strLengths[i], 0);
		m_klcpXY[i][1].resize(m_strLengths[i], m_kv);
		// last end
		for(int j = 1; j < m_kv; j++)
		{
			auto xit = m_strLengths[i] - j;
			if(xit >= 0)
			m_klcpXY[i][1][xit] = j;
		}
	}
	// get all the internal nodes
	std::vector<InternalNode> uNodes;
	selectInternalNodes0(uNodes);

	// for each internal node
	for(auto nit = uNodes.begin(); nit != uNodes.end(); nit++)
	{
		#ifdef DEBUG
		(*nit).dwriteln(m_aCfg.lfs);
		#endif
		//   collect tuples for each position (going left and right)
		//      (i, i', 0/1) i' = gisa[gsa[i] + d + 1]
		std::vector<L1Suffix> choppedSfxs;
		chopPrefix0(*nit, choppedSfxs);
		// update lcp using sorted tuples using a double pass
		computeK(*nit, choppedSfxs, m_kv - 1);
    	}
}

void ExactLCPk::compute()
{
	m_kv = m_aCfg.kv > 0 ? m_aCfg.kv : 0;
	if (m_kv == 0)	
	{
		compute0();
	} 
	else 
	{
	        computeK();
	}
	m_aCfg.lfs << "  [" << m_nPass << ", "  << m_passSizes << "]," << std::endl;
}

void ExactLCPk::computeTest(int k)
{
	m_kv= k;
	if (m_kv == 0)
	{
		compute0();
    	} 
	else 
	{
        	computeK();
    	}
}
