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

#ifndef EXACTLCPK_H
#define EXACTLCPK_H

#include "defs.hpp"
#include "AppConfig.hpp"
#include "ReadsDB.hpp"
#include "rmq_support_sparse_table.hpp"
#include <stack>

using namespace std;

struct InternalNode{
    int32_t m_leftBound;
    int32_t m_rightBound;
    int32_t m_stringDepth;
    int32_t m_delta;

    bool operator< (const InternalNode& other) const {
        return
            (m_leftBound < other.m_leftBound) ? true :
            ((m_leftBound == other.m_leftBound) ? (m_rightBound < other.m_rightBound) :
             false);
    }

    bool operator== (const InternalNode& other) const{
        return
            (m_leftBound == other.m_leftBound) &&
            (m_rightBound == other.m_rightBound) &&
            (m_stringDepth == other.m_stringDepth);
    }

    InternalNode(){
        m_leftBound = m_rightBound = -1;
        m_stringDepth = m_delta = 0;
    }

    void write(std::ostream& ots, const char *sepStr = "\t") const{
        ots << m_leftBound << sepStr << m_rightBound << sepStr
            << m_stringDepth << sepStr << m_delta;
    };

    void writeln(std::ostream& ots) const{
        write(ots);
        ots << std::endl;
    }

    void dwriteln(std::ostream& ots) const{
        ots << " [";
        write(ots, ",\t");
        ots << "]," << std::endl;
    }

};

struct L1Suffix{
  //      (c, c', 0/1), where c' = gisa[gsa[c] + d + 1]
    int32_t m_startPos; // starting position
    int32_t m_errSAPos; // position after one error's SA loc.
    int32_t m_srcStr;   // source string

    bool operator< (const L1Suffix& other) const {
        return( m_errSAPos < other.m_errSAPos );
    }

    L1Suffix(int32_t spos, int32_t epos, int32_t src):
        m_startPos(spos), m_errSAPos(epos), m_srcStr(src) { }

    L1Suffix(){}

    void write(std::ostream& ots, const char *sepStr = "\t") const{
        ots << m_startPos << sepStr << m_errSAPos << sepStr
            << m_srcStr;
    }

    void emit(std::ostream& ots) const{
        ots << m_startPos << "\t" << m_srcStr;
    }

    void writeln(std::ostream& ots) const{
        write(ots);
        ots << std::endl;
    }

    void dwrite(std::ostream& ots) const{
        ots << "\t[";
        write(ots, ",\t");
        ots << "],";
    }

    void dwriteln(std::ostream& ots) const{
        dwrite(ots);
        ots << std::endl;
    }

};

class UpperBoundCheck{
public:
    bool operator()(const int32_t& value, const int32_t& bound) const{
        return (value < bound);
    }
};

class LowerBoundCheck{
public:
    bool operator()(const int32_t& value, const int32_t& bound) const{
        return (value > bound);
    }
};

class IncrPointer{
public:
    void operator()(int32_t& value){
        value += 1;
    }
};

class DecrPointer{
public:
    void operator()(int32_t& value){
        value -= 1;
    }
};

class ExactLCPk{
private:
    AppConfig& m_aCfg;
    ivec_t m_klcpXY[2][2];
    ivec_t m_klcpHistoXY[2];
    int32_t m_strLengths[2];
    int32_t m_strLenPfx[2];
    int32_t m_shiftPos[2];
    std::string m_strXY;
    ivec_t m_gsa, m_gisa, m_glcp;
    int m_kv;
    double m_nPass;
    double m_passSizes;

    rmq_support_sparse_table<ivec_t, true, ivec_t> m_rangeMinQuery;

    int32_t leftBound0(int32_t curLeaf);
    int32_t rightBound0(int32_t curLeaf);
    void selectInternalNodes0(std::vector<InternalNode>& uNodes);
    void chopPrefix0(const InternalNode& uNode,
                     std::vector<L1Suffix>& leaves);

    void updateExactLCPk(InternalNode& uNode, std::vector<L1Suffix>& leaves);

    void eliminateDupes(std::vector<InternalNode>& uNodes);

    inline int32_t rangeMinLCP(const int32_t& t1, const int32_t& t2){
        if(t1 < 0 || t2 < 0)
            return 0;
        int32_t mxv = std::max(t1, t2);
        if(mxv >= (int32_t)m_gsa.size())
            return 0;
        int32_t mnv = std::min(t1, t2);
        int32_t rpos = m_rangeMinQuery(mnv == mxv ? mnv : (mnv + 1), mxv);
        return m_glcp[rpos];
    }

    inline int32_t rangeMinLCP(const L1Suffix& m1, const L1Suffix& m2){
        return rangeMinLCP(m1.m_errSAPos, m2.m_errSAPos);
    }

    inline int32_t sufRangeMinLCP(int32_t t1, int32_t t2) const{
        if(t1 < 0 || t2 < 0)
            return 0;
        if(std::max(t1, t2) >= (int32_t)m_gsa.size())
            return 0;
        int32_t st1 = m_gisa[t1];
        int32_t st2 = m_gisa[t2];
        int32_t mxv = std::max(st1, st2);
        int32_t mnv = std::min(st1, st2);
        int32_t rpos = m_rangeMinQuery(mnv == mxv ? mnv : (mnv + 1), mxv);
        return m_glcp[rpos];
    }

    inline int32_t updatePassLCP(const int32_t& t1, const int32_t& t2){
        if(t1 < 0 || t2 < 0)
            return 0;
        int32_t mxv = std::max(t1, t2);
        if(mxv > (int32_t)m_gsa.size())
            return 0;
        int32_t mnv = std::min(t1, t2);
        int32_t rpos = m_rangeMinQuery(mnv == mxv ? mnv : (mnv + 1), mxv);
        return 1 + m_glcp[rpos];
    }

    inline int32_t updatePassLCP(const L1Suffix& m1, const L1Suffix& m2){
        return updatePassLCP(m1.m_errSAPos, m2.m_errSAPos);
    }

    inline int32_t strPos(const InternalNode& uNode, const L1Suffix& sfx){
        return sfx.m_startPos - uNode.m_delta - m_shiftPos[sfx.m_srcStr];
    }

    template<typename BoundChecker, typename NextPointer>
    void updatePass(int32_t src_ptr, int32_t tgt_ptr,
                    const int32_t& tgt_bound,
                    const InternalNode& uNode,
                    const std::vector<L1Suffix>& leaves
#ifdef DEBUG
                    ,
                    const std::string& dbgStr
#endif
                    ) {
        BoundChecker bound_check;
        NextPointer next_ptr;
        // move the pointer until we reach the first src, target
        while(bound_check(tgt_ptr, tgt_bound)){
            if(leaves[tgt_ptr].m_srcStr == leaves[src_ptr].m_srcStr)
                break;
            src_ptr = tgt_ptr;
            next_ptr(tgt_ptr);
        }
        // if not within the bounds leave
        if(!bound_check(tgt_ptr, tgt_bound) )
            return;
        while(true){
            int32_t rmin = 0;
            int32_t tgt = leaves[tgt_ptr].m_srcStr,
                tpos = strPos(uNode, leaves[tgt_ptr]);
            // - get LCP between src_ptr and tgt_ptr from RMQ
            rmin = updatePassLCP(leaves[src_ptr], leaves[tgt_ptr]);
		
            int32_t score = uNode.m_stringDepth + uNode.m_delta + rmin;
#ifdef DEBUG
            m_aCfg.lfs << "\t[ \""
                       << dbgStr << "\",\t"
                       << score << ",\t"
                       << tgt << ",\t"
                       << tpos << ",\t"
                       << uNode.m_stringDepth << ",\t"
                       << uNode.m_delta << ",\t"
                       << rmin << ",\t"
                ;
            leaves[src_ptr].write(m_aCfg.lfs, ",\t");
            m_aCfg.lfs << ",\t";
            leaves[tgt_ptr].write(m_aCfg.lfs, ",\t");
            m_aCfg.lfs  << "]," << std::endl;
#endif
            assert(tpos >= 0);
            assert(tpos < (int32_t)m_klcpXY[tgt][1].size());
            // - update target's LCP, if score is higher
            if(score > m_klcpXY[tgt][1][tpos]){
                m_klcpXY[tgt][0][tpos] = strPos(uNode, leaves[src_ptr]);
                m_klcpXY[tgt][1][tpos] = score;
            }
            // - update tgt_ptr; quit if out of bounds
            int32_t prev_tgt = tgt_ptr;
            next_ptr(tgt_ptr);
            if(!bound_check(tgt_ptr, tgt_bound))
                break;
            // - update src_ptr, if tgt_ptr switches string source
            //if(leaves[tgt_ptr].m_srcStr == leaves[src_ptr].m_srcStr)
                src_ptr = prev_tgt;
        }
    }

    void computeK(InternalNode& uNode, std::vector<L1Suffix>& uLeaves,
                  int searchLevel);
    void computeK();
    int32_t leftBoundK(const std::vector<L1Suffix>& trieLeaves,
                       int32_t curLeaf);
    int32_t rightBoundK(const std::vector<L1Suffix>& trieLeaves,
                        int32_t curLeaf);
    void selectInternalNodesK(const InternalNode& prevNode,
                              const std::vector<L1Suffix>& leaves,
                              std::vector<InternalNode>& trieNodes);
    void chopPrefixK(const InternalNode& uNode,
                     const std::vector<L1Suffix>& uLeaves,
                     std::vector<L1Suffix>& trieLeaves);
    void compute0();
    void selectSuffixes0(const InternalNode& uNode,
                         std::vector<L1Suffix>& leaves);

    int32_t ansvLB0(int32_t curLeaf, int32_t ansvBound);
    int32_t ansvRB0(int32_t curLeaf, int32_t ansvBound);
    int32_t ansvLBK(const std::vector<L1Suffix>& trieLeaves,
                    int32_t curLeaf, int32_t minIdx,
                    int32_t ansvBound);
    int32_t ansvRBK(const std::vector<L1Suffix>& trieLeaves,
                    int32_t curLeaf, int32_t ansvBound);
public:
    friend class HeuristicLCPk;
    friend class LRHeuristicLCPk;
    ExactLCPk(const std::string& x,
              AppConfig& cfg);
    void print(std::ostream& ofs);
    void compute();
    auto getkLCP() -> const ivec_t (&)[2][2] {
        return m_klcpXY;
    }

    auto getkLCPHisto() -> const ivec_t (&)[2] {
      return m_klcpHistoXY;
    }
    void computeTest(int k);
};

template<typename SizeType, typename CountType>
void ansvBounds0(std::vector<SizeType>& sa,
                 std::vector<CountType>& lcp,
                 std::vector<SizeType>& left_bounds,
                 std::vector<SizeType>& right_bounds){
    std::stack<SizeType> s;
    // find right and left snv
    left_bounds.resize(sa.size());
    for(SizeType ix = 0; ix < (SizeType) sa.size(); ix++){
        while(!s.empty() &&
              lcp[s.top()] >= lcp[ix])
            s.pop();
        if(s.empty())
            left_bounds[ix] = 0;
        else
            left_bounds[ix] = s.top();
        s.push(ix);
    }

    while(!s.empty()) s.pop();

    right_bounds.resize(sa.size());
    for (SizeType ix = sa.size(); ix > 0; --ix) {
        while (!s.empty() &&
               lcp[s.top()] >= lcp[ix - 1])
            s.pop();
        if(s.empty())
            right_bounds[ix - 1] = sa.size() - 1;
        else
            right_bounds[ix - 1] = (s.top() > 0) ? s.top() - 1 : 0;
        s.push(ix - 1);
    }
}

template<typename SizeType, typename CountType,
         typename RangeMinLCP>
void ansvBoundsK(const std::vector<L1Suffix>& trieLeaves,
                 RangeMinLCP& rmin_qry,
                 std::vector<SizeType>& left_bounds,
                 std::vector<SizeType>& right_bounds){
    if(trieLeaves.size() == 0)
        return;
    std::stack<SizeType> s;
    // todo: find right and left snv
    left_bounds.resize(trieLeaves.size());
    left_bounds[0] = 0;
    for(SizeType ix = 1; ix < (SizeType)trieLeaves.size(); ix++){
        while(!s.empty() && // m_glcp[s.top()] >= m_glcp[ix]
              (rmin_qry(trieLeaves[s.top() - 1], trieLeaves[s.top()])
               >=
               rmin_qry(trieLeaves[ix - 1], trieLeaves[ix])))
            s.pop();
        if(s.empty())
            left_bounds[ix] = 0;
        else
            left_bounds[ix] = s.top();
        s.push(ix);
    }

    while(!s.empty()) s.pop();

    right_bounds.resize(trieLeaves.size());
    for (SizeType ix = trieLeaves.size(); ix > 1; --ix) {
        while (!s.empty() && // m_glcp[s.top()] >= m_glcp[ix - 1]
               (rmin_qry(trieLeaves[s.top() - 1], trieLeaves[s.top()])
                >=
                rmin_qry(trieLeaves[ix - 2], trieLeaves[ix - 1])))
            s.pop();
        if(s.empty())
            right_bounds[ix - 1] = trieLeaves.size() - 1;
        else
            right_bounds[ix - 1] = (s.top() > 0) ? s.top() - 1 : 0;
        s.push(ix - 1);
    }
}



#endif /* EXACTLCPK_H */
