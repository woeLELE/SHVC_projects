/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2015, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include "TEncTop.h"
#include "TEncCu.h"
#include "TEncAnalyze.h"
#include "TLibCommon/Debug.h"
#include "GMM.h"
#include "ITEW.h"
#include <cmath>
#include <algorithm>
using namespace std;

class CM
{
public:
	CM() {}
	CM(p c0_, vector<p> c1_, vector<p> c2_) : c0(c0_), c1(c1_), c2(c2_) {}
	bool isValid;
	p c0;
	vector<p> c1;
	vector<p> c2;

};

class CG
{
public:
	CG() {}
	CG(vector<vector<dg>> lv1_, vector<vector<dg>> lv2_) : lv1(lv1_), lv2(lv2_) {}
	bool isValid;
	vector<vector<dg>> lv1;
	vector<vector<dg>> lv2;
};

const double cond = 0.5;
const double fcond = 0.95;

vector<vector<dg>> tlv1;
vector<vector<dg>> tlv2;

vector<double> paras1({ 0.9976, -0.0825, -0.06541, 0.016, -0.001347, 3.888e-05 });
vector<double> paras2({ 0.9911, -0.002199, -0.01533, 0.001641, -5.38e-05 });

int strides[] = { 16, 8, 4 };
double RD_direct_ILR, RD_direct_merge, a = 0, b = 0, c = 0, m = 0;
const double limset = 0.95;
int m0 = 0;
int n0 = 0;
int ctuIndex = 0;
int fr = 0;
int px = 0, py = 0;
int px1 = 0, py1 = 0;
int px2 = 0, py2 = 0;
int px3 = 0, py3 = 0;
int indexD1 = 0;
int indexD2 = 0;
vector<int> pos;
bool gotSize = false;
bool inited = false;

typedef pair<int, double> p;
typedef pair<int, int> pp;
vector<vector<p>> curCosts0;
vector<vector<p>> curCosts1;
vector<vector<p>> curCosts2;
vector<vector<p>> curCosts3;
vector<vector<p>> preCosts0;
vector<vector<p>> preCosts1;
vector<vector<p>> preCosts2;
vector<vector<p>> preCosts3;

vector<vector<CM>> curmsg;
vector<vector<CM>> premsg;

vector<vector<CG>> curG;

vector<pp> prePos({ pp(0, -1), pp(0, 0), pp(-1, -1), pp(-1, 0), pp(-1, +1) });
vector<pp> curPos({ pp(0, -1), pp(-1, -1), pp(-1, 0), pp(-1, 1) });
vector<pp> curPos2({ pp(0, -1), pp(-1, 0) });

vector<int> index(4, 0);

//ofstream out0("输出0.txt");
//ofstream out1("输出1.txt");
//ofstream out2("输出2.txt");
//ofstream out3("输出3.txt");
//ofstream tout("调试输出.txt");
p cost0;
vector<p> cost1;
vector<p> cost2;
bool cease = false;
vector<double> pf(3, 0);
//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

/**
 \param    uhTotalDepth  total number of allowable depth
 \param    uiMaxWidth    largest CU width
 \param    uiMaxHeight   largest CU height
 \param    chromaFormat  chroma format
 */
Void TEncCu::create(UChar uhTotalDepth, UInt uiMaxWidth, UInt uiMaxHeight, ChromaFormat chromaFormat)
{
	Int i;

	m_uhTotalDepth = uhTotalDepth + 1;
	m_ppcBestCU = new TComDataCU * [m_uhTotalDepth - 1];
	m_ppcTempCU = new TComDataCU * [m_uhTotalDepth - 1];

	m_ppcPredYuvBest = new TComYuv * [m_uhTotalDepth - 1];
	m_ppcResiYuvBest = new TComYuv * [m_uhTotalDepth - 1];
	m_ppcRecoYuvBest = new TComYuv * [m_uhTotalDepth - 1];
	m_ppcPredYuvTemp = new TComYuv * [m_uhTotalDepth - 1];
	m_ppcResiYuvTemp = new TComYuv * [m_uhTotalDepth - 1];
	m_ppcRecoYuvTemp = new TComYuv * [m_uhTotalDepth - 1];
	m_ppcOrigYuv = new TComYuv * [m_uhTotalDepth - 1];

	UInt uiNumPartitions;
	for (i = 0; i < m_uhTotalDepth - 1; i++)
	{
		uiNumPartitions = 1 << ((m_uhTotalDepth - i - 1) << 1);
		UInt uiWidth = uiMaxWidth >> i;
		UInt uiHeight = uiMaxHeight >> i;

		m_ppcBestCU[i] = new TComDataCU;
		m_ppcBestCU[i]->create(chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1));
		m_ppcTempCU[i] = new TComDataCU;
		m_ppcTempCU[i]->create(chromaFormat, uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1));

		m_ppcPredYuvBest[i] = new TComYuv;
		m_ppcPredYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);
		m_ppcResiYuvBest[i] = new TComYuv;
		m_ppcResiYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);
		m_ppcRecoYuvBest[i] = new TComYuv;
		m_ppcRecoYuvBest[i]->create(uiWidth, uiHeight, chromaFormat);

		m_ppcPredYuvTemp[i] = new TComYuv;
		m_ppcPredYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);
		m_ppcResiYuvTemp[i] = new TComYuv;
		m_ppcResiYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);
		m_ppcRecoYuvTemp[i] = new TComYuv;
		m_ppcRecoYuvTemp[i]->create(uiWidth, uiHeight, chromaFormat);

		m_ppcOrigYuv[i] = new TComYuv;
		m_ppcOrigYuv[i]->create(uiWidth, uiHeight, chromaFormat);
	}

	m_bEncodeDQP = false;
	m_stillToCodeChromaQpOffsetFlag = false;
	m_cuChromaQpOffsetIdxPlus1 = 0;
	m_bFastDeltaQP = false;

	// initialize partition order.
	UInt* piTmp = &g_auiZscanToRaster[0];
	initZscanToRaster(m_uhTotalDepth, 1, 0, piTmp);
	initRasterToZscan(uiMaxWidth, uiMaxHeight, m_uhTotalDepth);

	// initialize conversion matrix from partition index to pel
	initRasterToPelXY(uiMaxWidth, uiMaxHeight, m_uhTotalDepth);
	f = fopen("1.txt", "a+");
}

Void TEncCu::destroy()
{
	Int i;

	for (i = 0; i < m_uhTotalDepth - 1; i++)
	{
		if (m_ppcBestCU[i])
		{
			m_ppcBestCU[i]->destroy();
			delete m_ppcBestCU[i];
			m_ppcBestCU[i] = NULL;
		}
		if (m_ppcTempCU[i])
		{
			m_ppcTempCU[i]->destroy();
			delete m_ppcTempCU[i];
			m_ppcTempCU[i] = NULL;
		}
		if (m_ppcPredYuvBest[i])
		{
			m_ppcPredYuvBest[i]->destroy();
			delete m_ppcPredYuvBest[i];
			m_ppcPredYuvBest[i] = NULL;
		}
		if (m_ppcResiYuvBest[i])
		{
			m_ppcResiYuvBest[i]->destroy();
			delete m_ppcResiYuvBest[i];
			m_ppcResiYuvBest[i] = NULL;
		}
		if (m_ppcRecoYuvBest[i])
		{
			m_ppcRecoYuvBest[i]->destroy();
			delete m_ppcRecoYuvBest[i];
			m_ppcRecoYuvBest[i] = NULL;
		}
		if (m_ppcPredYuvTemp[i])
		{
			m_ppcPredYuvTemp[i]->destroy();
			delete m_ppcPredYuvTemp[i];
			m_ppcPredYuvTemp[i] = NULL;
		}
		if (m_ppcResiYuvTemp[i])
		{
			m_ppcResiYuvTemp[i]->destroy();
			delete m_ppcResiYuvTemp[i];
			m_ppcResiYuvTemp[i] = NULL;
		}
		if (m_ppcRecoYuvTemp[i])
		{
			m_ppcRecoYuvTemp[i]->destroy();
			delete m_ppcRecoYuvTemp[i];
			m_ppcRecoYuvTemp[i] = NULL;
		}
		if (m_ppcOrigYuv[i])
		{
			m_ppcOrigYuv[i]->destroy();
			delete m_ppcOrigYuv[i];
			m_ppcOrigYuv[i] = NULL;
		}
	}
	if (m_ppcBestCU)
	{
		delete[] m_ppcBestCU;
		m_ppcBestCU = NULL;
	}
	if (m_ppcTempCU)
	{
		delete[] m_ppcTempCU;
		m_ppcTempCU = NULL;
	}

	if (m_ppcPredYuvBest)
	{
		delete[] m_ppcPredYuvBest;
		m_ppcPredYuvBest = NULL;
	}
	if (m_ppcResiYuvBest)
	{
		delete[] m_ppcResiYuvBest;
		m_ppcResiYuvBest = NULL;
	}
	if (m_ppcRecoYuvBest)
	{
		delete[] m_ppcRecoYuvBest;
		m_ppcRecoYuvBest = NULL;
	}
	if (m_ppcPredYuvTemp)
	{
		delete[] m_ppcPredYuvTemp;
		m_ppcPredYuvTemp = NULL;
	}
	if (m_ppcResiYuvTemp)
	{
		delete[] m_ppcResiYuvTemp;
		m_ppcResiYuvTemp = NULL;
	}
	if (m_ppcRecoYuvTemp)
	{
		delete[] m_ppcRecoYuvTemp;
		m_ppcRecoYuvTemp = NULL;
	}
	if (m_ppcOrigYuv)
	{
		delete[] m_ppcOrigYuv;
		m_ppcOrigYuv = NULL;
	}
	fclose(f);
}

/** \param    pcEncTop      pointer of encoder class
 */
Void TEncCu::init(TEncTop* pcEncTop)
{
	m_pcEncCfg = pcEncTop;
	m_pcPredSearch = pcEncTop->getPredSearch();
	m_pcTrQuant = pcEncTop->getTrQuant();
	m_pcRdCost = pcEncTop->getRdCost();

#if SVC_EXTENSION
	m_ppcTEncTop = pcEncTop->getLayerEnc();
#endif

	m_pcEntropyCoder = pcEncTop->getEntropyCoder();
	m_pcBinCABAC = pcEncTop->getBinCABAC();

	m_pppcRDSbacCoder = pcEncTop->getRDSbacCoder();
	m_pcRDGoOnSbacCoder = pcEncTop->getRDGoOnSbacCoder();

	m_pcRateCtrl = pcEncTop->getRateCtrl();
}

// ====================================================================================================================
// Public member functions
// Public member functions
// ====================================================================================================================

/**
 \param  pCtu pointer of CU data class
 */
Void TEncCu::compressCtu(TComDataCU* pCtu)
{
	// initialize CU data
	m_ppcBestCU[0]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
	m_ppcTempCU[0]->initCtu(pCtu->getPic(), pCtu->getCtuRsAddr());
	//index1 = 0;
	//index2 = 0;
#if N0383_IL_CONSTRAINED_TILE_SETS_SEI
	m_disableILP = xCheckTileSetConstraint(pCtu);
	m_pcPredSearch->setDisableILP(m_disableILP);
#endif
	if (m_ppcBestCU[0]->getPic()->getLayerId() > 0)
	{
		index.assign(4, 0);

		cost0 = p(-1, 0);
		cost1.assign(4, p(-1, 0));
		cost2.assign(16, p(-1, 0));

		tlv1.assign(2, vector<dg>(2));
		tlv2.assign(4, vector<dg>(4));

		pf.assign(3, -1);
	}
	if (m_ppcBestCU[0]->getPic()->getLayerId() > 0 && !gotSize)
	{
		m0 = m_ppcBestCU[0]->getPic()->getFrameWidthInCtus();
		n0 = m_ppcBestCU[0]->getPic()->getFrameHeightInCtus();
		gotSize = true;
	}

	if (!inited && gotSize)
	{
		cout << "内部表优化" << endl;
		cout << cond << endl;
		cout << fcond << endl;
		pos.push_back(m0 * n0 - 1);
		pos.push_back(m0 * n0 - m0 - 1);
		pos.push_back(m0 * n0 - m0);
		pos.push_back(m0 * n0 - m0 + 1);
		curCosts0.assign(n0, vector<p>(m0, p(-1, -1)));
		curCosts1.assign(n0 * 2, vector<p>(m0 * 2, p(-1, -1)));
		curCosts2.assign(n0 * 4, vector<p>(m0 * 4, p(-1, -1)));
		curCosts3.assign(n0 * 8, vector<p>(m0 * 8, p(-1, -1)));

		curmsg.assign(n0, vector<CM>(m0, CM()));
		curG.assign(n0, vector<CG>(m0, CG()));
		//curmsg.assign(n0, vector<CM>(m0, CM(p(-1, 0), vector<p>(4, p(-1, 0)), vector<p>(16, p(-1, 0)))));
		inited = true;
	}


	if (m0 > 0 && n0 > 0 && m_ppcBestCU[0]->getPic()->getLayerId() > 0)
	{
		px = ctuIndex / m0;
		py = ctuIndex % m0;
	}

	// analysis of CU
	DEBUG_STRING_NEW(sDebug)

		xCompressCU(m_ppcBestCU[0], m_ppcTempCU[0], 0 DEBUG_STRING_PASS_INTO(sDebug));
	DEBUG_STRING_OUTPUT(std::cout, sDebug)
		int max = 0, i, j, x, y, d[16][16], m = -1, n, k = 0;
#if ADAPTIVE_QP_SELECTION
	if (m_pcEncCfg->getUseAdaptQpSelect())
	{
		if (pCtu->getSlice()->getSliceType() != I_SLICE) //IIII
		{
			xCtuCollectARLStats(pCtu);
		}
	}
#endif
	if (m_ppcBestCU[0]->getPic()->getLayerId() > 0)
	{
		if (index[1] == 4 && index[2] == 16)
		{
			int iCount = 0;
			int iWidthInPart = MAX_CU_SIZE >> 2;
			for (int i = 0; i < pCtu->getTotalNumPart(); i++)
			{
				if ((iCount & (iWidthInPart - 1)) == 0)
				{
					m++;
					n = 0;
				}
				d[m][n++] = pCtu->getDepth(g_auiRasterToZscan[i]);
				iCount++;
			}
			for (int depth = 0; depth < 3; depth++)		// 得到不同大小CU的深度
			{
				int stride = strides[depth];
				int idx = 0;
				for (i = 0; i < 16; i = i + stride)
				{
					for (j = 0; j < 16; j = j + stride)
					{
						max = 0;
						for (m = i; m < i + stride; m++)
						{
							for (n = j; n < j + stride; n++)
							{
								if (max < d[m][n])
									max = d[m][n];
								if (max > depth)
									break;
							}
							if (max > depth)
								break;
						}
						if (depth == 0)
							cost0.first = max > depth ? 1 : 0;
						else if (depth == 1)
							cost1[idx++].first = max > depth ? 1 : 0;
						else if (depth == 2)
							cost2[idx++].first = max > depth ? 1 : 0;
					}
				}
			}
			curmsg[px][py] = CM(cost0, cost1, cost2);
			curmsg[px][py].isValid = true;
			curG[px][py] = CG(tlv1, tlv2);
			curG[px][py].isValid = true;
		}
		else
		{
			curmsg[px][py].isValid = false;
			curG[px][py].isValid = false;
		}

		ctuIndex++;
	}
	if (ctuIndex == m0 * n0 && m_ppcBestCU[0]->getPic()->getLayerId() > 0)
	{
		ctuIndex = 0;
		fr++;
		preCosts0 = curCosts0;
		preCosts1 = curCosts1;
		preCosts2 = curCosts2;
		preCosts3 = curCosts3;

		premsg = curmsg;

		curCosts0.assign(n0, vector<p>(m0, p(-1, -1)));
		curCosts1.assign(n0 * 2, vector<p>(m0 * 2, p(-1, -1)));
		curCosts2.assign(n0 * 4, vector<p>(m0 * 4, p(-1, -1)));
		curCosts3.assign(n0 * 8, vector<p>(m0 * 8, p(-1, -1)));

		//curmsg.assign(n0, vector<CM>(m0, CM(p(-1, 0), vector<p>(4, p(-1, 0)), vector<p>(16, p(-1, 0)))));
		curmsg.assign(n0, vector<CM>(m0, CM()));
		curG.assign(n0, vector<CG>(m0, CG()));
	}
#if N0383_IL_CONSTRAINED_TILE_SETS_SEI
	xVerifyTileSetConstraint(pCtu);
#endif
}
/** \param  pCtu  pointer of CU data class
 */
Void TEncCu::encodeCtu(TComDataCU* pCtu)
{
	if (pCtu->getSlice()->getPPS()->getUseDQP())
	{
		setdQPFlag(true);
	}

	if (pCtu->getSlice()->getUseChromaQpAdj())
	{
		setCodeChromaQpAdjFlag(true);
	}

	// Encode CU data
	xEncodeCU(pCtu, 0, 0);
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================
//! Derive small set of test modes for AMP encoder speed-up
#if AMP_ENC_SPEEDUP
#if AMP_MRG
Void TEncCu::deriveTestModeAMP(TComDataCU* pcBestCU, PartSize eParentPartSize, Bool& bTestAMP_Hor, Bool& bTestAMP_Ver, Bool& bTestMergeAMP_Hor, Bool& bTestMergeAMP_Ver)
#else
Void TEncCu::deriveTestModeAMP(TComDataCU* pcBestCU, PartSize eParentPartSize, Bool& bTestAMP_Hor, Bool& bTestAMP_Ver)
#endif
{
	if (pcBestCU->getPartitionSize(0) == SIZE_2NxN)
	{
		bTestAMP_Hor = true;
	}
	else if (pcBestCU->getPartitionSize(0) == SIZE_Nx2N)
	{
		bTestAMP_Ver = true;
	}
	else if (pcBestCU->getPartitionSize(0) == SIZE_2Nx2N && pcBestCU->getMergeFlag(0) == false && pcBestCU->isSkipped(0) == false)
	{
		bTestAMP_Hor = true;
		bTestAMP_Ver = true;
	}

#if AMP_MRG
	//! Utilizing the partition size of parent PU
	if (eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N)
	{
		bTestMergeAMP_Hor = true;
		bTestMergeAMP_Ver = true;
	}

	if (eParentPartSize == NUMBER_OF_PART_SIZES) //! if parent is intra
	{
		if (pcBestCU->getPartitionSize(0) == SIZE_2NxN)
		{
			bTestMergeAMP_Hor = true;
		}
		else if (pcBestCU->getPartitionSize(0) == SIZE_Nx2N)
		{
			bTestMergeAMP_Ver = true;
		}
	}

	if (pcBestCU->getPartitionSize(0) == SIZE_2Nx2N && pcBestCU->isSkipped(0) == false)
	{
		bTestMergeAMP_Hor = true;
		bTestMergeAMP_Ver = true;
	}

	if (pcBestCU->getWidth(0) == 64)
	{
		bTestAMP_Hor = false;
		bTestAMP_Ver = false;
	}
#else
	//! Utilizing the partition size of parent PU
	if (eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N)
	{
		bTestAMP_Hor = true;
		bTestAMP_Ver = true;
	}

	if (eParentPartSize == SIZE_2Nx2N)
	{
		bTestAMP_Hor = false;
		bTestAMP_Ver = false;
	}
#endif
}
#endif

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================
/** Compress a CU block recursively with enabling sub-CTU-level delta QP
 *  - for loop of QP value to compress the current CU with all possible QP
*/
#if AMP_ENC_SPEEDUP
Void TEncCu::xCompressCU(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth DEBUG_STRING_FN_DECLARE(sDebug_), PartSize eParentPartSize)
#else
Void TEncCu::xCompressCU(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, const UInt uiDepth)
#endif
{
	cease = false;
	TComPic* pcPic = rpcBestCU->getPic();
	//if (pcPic->getLayerId() > 0)
	//{
	//    out2 << pcPic->getFrameWidthInCtus() << " : " << pcPic->getFrameHeightInCtus() << endl;
	//}
	DEBUG_STRING_NEW(sDebug)
		const TComPPS& pps = *(rpcTempCU->getSlice()->getPPS());
	const TComSPS& sps = *(rpcTempCU->getSlice()->getSPS());

	// These are only used if getFastDeltaQp() is true
	const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight() >> sps.getLog2DiffMaxMinCodingBlockSize(), sps.getMaxCUHeight(), 32u);

	// get Original YUV data from picture
	m_ppcOrigYuv[uiDepth]->copyFromPicYuv(pcPic->getPicYuvOrg(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu());

	// variable for Cbf fast mode PU decision
	Bool doNotBlockPu = true;
	Bool earlyDetectionSkipMode = false;

	const UInt uiLPelX = rpcBestCU->getCUPelX();
	const UInt uiRPelX = uiLPelX + rpcBestCU->getWidth(0) - 1;
	const UInt uiTPelY = rpcBestCU->getCUPelY();
	const UInt uiBPelY = uiTPelY + rpcBestCU->getHeight(0) - 1;
	const UInt uiWidth = rpcBestCU->getWidth(0);

	Int iBaseQP = xComputeQP(rpcBestCU, uiDepth);
	Int iMinQP;
	Int iMaxQP;
	Bool isAddLowestQP = false;

	const UInt numberValidComponents = rpcBestCU->getPic()->getNumberValidComponents();

	dg cc;
	if ((uiDepth == 1 || uiDepth == 2) && pcPic->getLayerId() > 0)
	{
		int i, j;
		float s = 0;
		UInt uiPartSize = rpcBestCU->getWidth(0);
		const ComponentID compID = ComponentID(0);
		const Int Width = uiPartSize >> m_ppcOrigYuv[uiDepth]->getComponentScaleX(compID);   // 整个CTU的边长
		const Pel* pSrc0 = m_ppcOrigYuv[uiDepth]->getAddr(compID, 0, Width);
		const Int  iSrc0Stride = m_ppcOrigYuv[uiDepth]->getStride(compID);
		vector<vector<int>> p(Width, vector<int>(Width));
		double avg = 0;
		// 得到灰度图
		for (i = 0; i < Width; i++)
		{
			for (j = 0; j < Width; j++)
			{
				p[i][j] = pSrc0[j];
				avg += p[i][j];
				s += pSrc0[j];
			}
			pSrc0 += iSrc0Stride;
		}
		avg /= (double)(Width * Width);
		cc.first = avg;
		cc.second = p;
		if (isEdge(uiDepth, index[uiDepth]))
		{
			if (uiDepth == 1)
			{
				tlv1[index[1] / 2][index[1] % 2].first = avg;
				tlv1[index[1] / 2][index[1] % 2].second = p;
			}
			else if (uiDepth == 2)
			{
				tlv2[index[2] / 4][index[2] % 4].first = avg;
				tlv2[index[2] / 4][index[2] % 4].second = p;
			}
		}
	}
	if (uiDepth <= pps.getMaxCuDQPDepth())
	{
		Int idQP = m_pcEncCfg->getMaxDeltaQP();
		iMinQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP - idQP);
		iMaxQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP + idQP);
	}
	else
	{
		iMinQP = rpcTempCU->getQP(0);
		iMaxQP = rpcTempCU->getQP(0);
	}

	if (m_pcEncCfg->getUseRateCtrl())
	{
		iMinQP = m_pcRateCtrl->getRCQP();
		iMaxQP = m_pcRateCtrl->getRCQP();
	}

	// transquant-bypass (TQB) processing loop variable initialisation ---

	const Int lowestQP = iMinQP; // For TQB, use this QP which is the lowest non TQB QP tested (rather than QP'=0) - that way delta QPs are smaller, and TQB can be tested at all CU levels.

	if ((pps.getTransquantBypassEnableFlag()))
	{
		isAddLowestQP = true; // mark that the first iteration is to cost TQB mode.
		iMinQP = iMinQP - 1;  // increase loop variable range by 1, to allow testing of TQB mode along with other QPs
		if (m_pcEncCfg->getCUTransquantBypassFlagForceValue())
		{
			iMaxQP = iMinQP;
		}
	}

	TComSlice* pcSlice = rpcTempCU->getPic()->getSlice(rpcTempCU->getPic()->getCurrSliceIdx());
	const Bool bBoundary = !(uiRPelX < sps.getPicWidthInLumaSamples() && uiBPelY < sps.getPicHeightInLumaSamples());

	if (!bBoundary)
	{
#if HIGHER_LAYER_IRAP_SKIP_FLAG
		if (m_pcEncCfg->getSkipPictureAtArcSwitch() && m_pcEncCfg->getAdaptiveResolutionChange() > 0 && pcSlice->getLayerId() == 1 && pcSlice->getPOC() == m_pcEncCfg->getAdaptiveResolutionChange())
		{
			Int iQP = iBaseQP;
			const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP);

			if (bIsLosslessMode)
			{
				iQP = lowestQP;
			}

			rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

			xCheckRDCostMerge2Nx2N(rpcBestCU, rpcTempCU, &earlyDetectionSkipMode, true);
		}
		else
		{
#endif
#if ENCODER_FAST_MODE
			Bool testInter = true;
			if (rpcBestCU->getPic()->getLayerId() > 0)
			{
				if (pcSlice->getSliceType() == P_SLICE && pcSlice->getNumRefIdx(REF_PIC_LIST_0) == pcSlice->getActiveNumILRRefIdx())
				{
					testInter = false;
				}
				if (pcSlice->getSliceType() == B_SLICE && pcSlice->getNumRefIdx(REF_PIC_LIST_0) == pcSlice->getActiveNumILRRefIdx() && pcSlice->getNumRefIdx(REF_PIC_LIST_1) == pcSlice->getActiveNumILRRefIdx())
				{
					testInter = false;
				}
			}
#endif
			for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
			{
				const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP);

				if (bIsLosslessMode)
				{
					iQP = lowestQP;
				}

				m_cuChromaQpOffsetIdxPlus1 = 0;
				if (pcSlice->getUseChromaQpAdj())
				{
					/* Pre-estimation of chroma QP based on input block activity may be performed
		 * here, using for example m_ppcOrigYuv[uiDepth] */
		 /* To exercise the current code, the index used for adjustment is based on
* block position
*/
					Int lgMinCuSize = sps.getLog2MinCodingBlockSize() +
						std::max<Int>(0, sps.getLog2DiffMaxMinCodingBlockSize() - Int(pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth()));
					m_cuChromaQpOffsetIdxPlus1 = ((uiLPelX >> lgMinCuSize) + (uiTPelY >> lgMinCuSize)) % (pps.getPpsRangeExtension().getChromaQpOffsetListLen() + 1);
				}

				rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

				// do inter modes, SKIP and 2Nx2N
				if (pcPic->getLayerId() > 0)
				{
#if ENCODER_FAST_MODE == 1
					if (rpcBestCU->getSlice()->getSliceType() != I_SLICE && testInter)
#else
					if (rpcBestCU->getSlice()->getSliceType() != I_SLICE)
#endif
					{
						for (Int refLayer = 0; refLayer < pcSlice->getActiveNumILRRefIdx(); refLayer++)
						{// 层间
							xCheckRDCostILRUni(rpcBestCU, rpcTempCU, pcSlice->getVPS()->getRefLayerId(pcSlice->getLayerId(), pcSlice->getInterLayerPredLayerIdc(refLayer)));
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
						}
						RD_direct_ILR = rpcBestCU->getTotalCost();

						// SKIP
						xCheckRDCostMerge2Nx2N(rpcBestCU, rpcTempCU DEBUG_STRING_PASS_INTO(sDebug), &earlyDetectionSkipMode); //by Merge for inter_2Nx2N
						rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
						RD_direct_merge = rpcBestCU->getTotalCost();


						// 2Nx2N
						if (m_pcEncCfg->getUseEarlySkipDetection())
						{
							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode); //by Competition for inter_2Nx2N
						}

#if ENCODER_FAST_MODE == 2
						if (testInter)
						{
#endif
							if (!m_pcEncCfg->getUseEarlySkipDetection())
							{
								// 2Nx2N, NxN
								xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
								rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
								if (m_pcEncCfg->getUseCbfFastMode())
								{
									doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
								}
							}
#if ENCODER_FAST_MODE == 2
						}
#endif
					}
				}
				else
				{
#if ENCODER_FAST_MODE == 1
					if (rpcBestCU->getSlice()->getSliceType() != I_SLICE && testInter)
#else
					if (rpcBestCU->getSlice()->getSliceType() != I_SLICE)
#endif
					{
						// 2Nx2N
						if (m_pcEncCfg->getUseEarlySkipDetection())
						{
							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode); //by Competition for inter_2Nx2N
						}
						// SKIP
						xCheckRDCostMerge2Nx2N(rpcBestCU, rpcTempCU DEBUG_STRING_PASS_INTO(sDebug), &earlyDetectionSkipMode); //by Merge for inter_2Nx2N
						rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

#if ENCODER_FAST_MODE == 2
						if (testInter)
						{
#endif
							if (!m_pcEncCfg->getUseEarlySkipDetection())
							{
								// 2Nx2N, NxN
								xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
								rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
								if (m_pcEncCfg->getUseCbfFastMode())
								{
									doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
								}
							}
#if ENCODER_FAST_MODE == 2
						}
#endif
					}
				}

				if (bIsLosslessMode) // Restore loop variable if lossless mode was searched.
				{
					iQP = iMinQP;
				}
			}

			if (!earlyDetectionSkipMode)
			{
				for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
				{
					const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP); // If lossless, then iQP is irrelevant for subsequent modules.

					if (bIsLosslessMode)
					{
						iQP = lowestQP;
					}

					rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

					// do inter modes, NxN, 2NxN, and Nx2N
#if ENCODER_FAST_MODE
					if (rpcBestCU->getSlice()->getSliceType() != I_SLICE && testInter)
#else
					if (rpcBestCU->getSlice()->getSliceType() != I_SLICE)
#endif
					{
						// 2Nx2N, NxN

						if (!((rpcBestCU->getWidth(0) == 8) && (rpcBestCU->getHeight(0) == 8)))
						{
							if (uiDepth == sps.getLog2DiffMaxMinCodingBlockSize() && doNotBlockPu)
							{
								xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug));
								rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
							}
						}

						if (doNotBlockPu)
						{
							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_Nx2N DEBUG_STRING_PASS_INTO(sDebug));
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
							if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_Nx2N)
							{
								doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
							}
						}
						if (doNotBlockPu)
						{
							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxN DEBUG_STRING_PASS_INTO(sDebug));
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
							if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxN)
							{
								doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
							}
						}

						//! Try AMP (SIZE_2NxnU, SIZE_2NxnD, SIZE_nLx2N, SIZE_nRx2N)
						if (sps.getUseAMP() && uiDepth < sps.getLog2DiffMaxMinCodingBlockSize())
						{
#if AMP_ENC_SPEEDUP
							Bool bTestAMP_Hor = false, bTestAMP_Ver = false;

#if AMP_MRG
							Bool bTestMergeAMP_Hor = false, bTestMergeAMP_Ver = false;

							deriveTestModeAMP(rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver, bTestMergeAMP_Hor, bTestMergeAMP_Ver);
#else
							deriveTestModeAMP(rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver);
#endif

							//! Do horizontal AMP
							if (bTestAMP_Hor)
							{
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnU DEBUG_STRING_PASS_INTO(sDebug));
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
									if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU)
									{
										doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
									}
								}
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnD DEBUG_STRING_PASS_INTO(sDebug));
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
									if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD)
									{
										doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
									}
								}
							}
#if AMP_MRG
							else if (bTestMergeAMP_Hor)
							{
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnU DEBUG_STRING_PASS_INTO(sDebug), true);
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
									if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU)
									{
										doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
									}
								}
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnD DEBUG_STRING_PASS_INTO(sDebug), true);
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
									if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD)
									{
										doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
									}
								}
							}
#endif

							//! Do horizontal AMP
							if (bTestAMP_Ver)
							{
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nLx2N DEBUG_STRING_PASS_INTO(sDebug));
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
									if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N)
									{
										doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
									}
								}
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nRx2N DEBUG_STRING_PASS_INTO(sDebug));
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
								}
							}
#if AMP_MRG
							else if (bTestMergeAMP_Ver)
							{
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nLx2N DEBUG_STRING_PASS_INTO(sDebug), true);
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
									if (m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N)
									{
										doNotBlockPu = rpcBestCU->getQtRootCbf(0) != 0;
									}
								}
								if (doNotBlockPu)
								{
									xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nRx2N DEBUG_STRING_PASS_INTO(sDebug), true);
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
								}
							}
#endif

#else
							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnU);
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_2NxnD);
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nLx2N);
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

							xCheckRDCostInter(rpcBestCU, rpcTempCU, SIZE_nRx2N);
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

#endif
						}
					}


					// test PCM
					if (sps.getUsePCM() && rpcTempCU->getWidth(0) <= (1 << sps.getPCMLog2MaxSize()) && rpcTempCU->getWidth(0) >= (1 << sps.getPCMLog2MinSize()))
					{
						UInt uiRawBits = getTotalBits(rpcBestCU->getWidth(0), rpcBestCU->getHeight(0), rpcBestCU->getPic()->getChromaFormat(), sps.getBitDepths().recon);
						UInt uiBestBits = rpcBestCU->getTotalBits();
						if ((uiBestBits > uiRawBits) || (rpcBestCU->getTotalCost() > m_pcRdCost->calcRdCost(uiRawBits, 0)))
						{
							xCheckIntraPCM(rpcBestCU, rpcTempCU);
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
						}
					}
#if ENCODER_FAST_MODE
#if N0383_IL_CONSTRAINED_TILE_SETS_SEI
					if (pcPic->getLayerId() > 0/* && !m_disableILP*/)
#else
					if (pcPic->getLayerId() > 0)
#endif
					{
						double res = 0.499;
						if (uiDepth == 1)
						{
							px1 = px * 2 + index[uiDepth] / 2;
							py1 = py * 2 + index[uiDepth] % 2;
						}
						else if (uiDepth == 2)
						{
							px2 = px * 4 + index[uiDepth] / 4;
							py2 = py * 4 + index[uiDepth] % 4;
						}
						else if (uiDepth == 3)
						{
							px3 = px * 8 + index[uiDepth] / 8;
							py3 = py * 8 + index[uiDepth] % 8;
						}

						Double ILRCost = 0, intraCost = 0.0;
						for (Int refLayer = 0; refLayer < pcSlice->getActiveNumILRRefIdx(); refLayer++)
						{// 层间
							xCheckRDCostILRUni(rpcBestCU, rpcTempCU, pcSlice->getVPS()->getRefLayerId(pcSlice->getLayerId(), pcSlice->getInterLayerPredLayerIdc(refLayer)));
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
						}
						// 在这里进行概率计算
						ILRCost = rpcBestCU->getTotalCost();
						if (fr > 0 && px > 0 && px < n0 - 1 && py > 0 && py < m0 - 1)
						{
							if (uiDepth == 0)
							{
								vector<double> rec0;
								vector<double> rec1;
								int x, y;
								double e0 = 0.0;
								double e1 = 0.0;
								for (auto e : prePos)
								{
									x = px + e.first;
									y = py + e.second;
									if (preCosts0[x][y].first == 0)
									{
										rec0.push_back(preCosts0[x][y].second);
										e0 += rec0.back();
									}
									else if (preCosts0[x][y].first == 1)
									{
										rec1.push_back(preCosts0[x][y].second);
										e1 += rec1.back();
									}
								}
								for (auto e : curPos)
								{
									x = px + e.first;
									y = py + e.second;
									if (curCosts0[x][y].first == 0)
									{
										rec0.push_back(curCosts0[x][y].second);
										e0 += rec0.back();
									}
									else if (curCosts0[x][y].first == 1)
									{
										rec1.push_back(curCosts0[x][y].second);
										e1 += rec1.back();
									}
								}
								if (rec0.size() < 1)
									res = 0;
								else if (rec1.size() < 1)
									res = 1;
								else
								{
									e0 /= rec0.size();
									e1 /= rec1.size();
									if (rec0.size() > 1 && rec1.size() > 1)
									{
										GMM gmm(rec0, rec1, e0, e1, ILRCost);
										gmm.init();
										gmm.train();
										res = gmm.getResult();
									}
								}
							}

							if (uiDepth == 1)
							{
								vector<double> rec0;
								vector<double> rec1;
								int x, y;
								double e0 = 0.0;
								double e1 = 0.0;
								for (auto e : prePos)
								{
									x = px1 + e.first;
									y = py1 + e.second;
									if (preCosts1[x][y].first == 0)
									{
										rec0.push_back(preCosts1[x][y].second);
										e0 += rec0.back();
									}
									else if (preCosts1[x][y].first == 1)
									{
										rec1.push_back(preCosts1[x][y].second);
										e1 += rec1.back();
									}
								}
								for (auto e : curPos)
								{
									x = px1 + e.first;
									y = py1 + e.second;
									if (curCosts1[x][y].first == 0)
									{
										rec0.push_back(curCosts1[x][y].second);
										e0 += rec0.back();
									}
									else if (curCosts1[x][y].first == 1)
									{
										rec1.push_back(curCosts1[x][y].second);
										e1 += rec1.back();
									}
								}

								if (rec0.size() < 1)
									res = 0;
								else if (rec1.size() < 1)
									res = 1;
								else
								{
									e0 /= rec0.size();
									e1 /= rec1.size();
									if (rec0.size() > 1 && rec1.size() > 1)
									{
										GMM gmm(rec0, rec1, e0, e1, ILRCost);
										gmm.init();
										gmm.train();
										res = gmm.getResult();
									}
								}
							}

							if (uiDepth == 2)
							{
								vector<double> rec0;
								vector<double> rec1;
								int x, y;
								double e0 = 0.0;
								double e1 = 0.0;
								for (auto e : prePos)
								{
									x = px2 + e.first;
									y = py2 + e.second;
									if (preCosts2[x][y].first == 0)
									{
										rec0.push_back(preCosts2[x][y].second);
										e0 += rec0.back();
									}
									else if (preCosts2[x][y].first == 1)
									{
										rec1.push_back(preCosts2[x][y].second);
										e1 += rec1.back();
									}
								}
								for (auto e : curPos)
								{
									x = px2 + e.first;
									y = py2 + e.second;
									if (curCosts2[x][y].first == 0)
									{
										rec0.push_back(curCosts2[x][y].second);
										e0 += rec0.back();
									}
									else if (curCosts2[x][y].first == 1)
									{
										rec1.push_back(curCosts2[x][y].second);
										e1 += rec1.back();
									}
								}
								if (rec0.size() < 1)
									res = 0;
								else if (rec1.size() < 1)
									res = 1;
								else
								{
									e0 /= rec0.size();
									e1 /= rec1.size();
									if (rec0.size() > 1 && rec1.size() > 1)
									{
										GMM gmm(rec0, rec1, e0, e1, ILRCost);
										gmm.init();
										gmm.train();
										res = gmm.getResult();
									}
								}
							}

							if (uiDepth == 3)
							{
								vector<double> rec0;
								vector<double> rec1;
								int x, y;
								double e0 = 0.0;
								double e1 = 0.0;
								for (auto e : prePos)
								{
									x = px3 + e.first;
									y = py3 + e.second;
									if (preCosts3[x][y].first == 0)
									{
										rec0.push_back(preCosts3[x][y].second);
										e0 += rec0.back();
									}
									else if (preCosts3[x][y].first == 1)
									{
										rec1.push_back(preCosts3[x][y].second);
										e1 += rec1.back();
									}
								}
								for (auto e : curPos)
								{
									x = px3 + e.first;
									y = py3 + e.second;
									if (curCosts3[x][y].first == 0)
									{
										rec0.push_back(curCosts3[x][y].second);
										e0 += rec0.back();
									}
									else if (curCosts3[x][y].first == 1)
									{
										rec1.push_back(curCosts3[x][y].second);
										e1 += rec1.back();
									}
								}
								if (rec0.size() < 1)
									res = 0;
								else if (rec1.size() < 1)
									res = 1;
								else
								{
									e0 /= rec0.size();
									e1 /= rec1.size();
									if (rec0.size() > 1 && rec1.size() > 1)
									{
										GMM gmm(rec0, rec1, e0, e1, ILRCost);
										gmm.init();
										gmm.train();
										res = gmm.getResult();
									}
								}
							}
						}

						if (res <= limset)
						{
							//Double intraCost = 0.0;

							if ((rpcBestCU->getSlice()->getSliceType() == I_SLICE) ||
#if ENCODER_FAST_MODE
								rpcBestCU->getPredictionMode(0) == NUMBER_OF_PREDICTION_MODES || // if there is no valid inter prediction
								!testInter ||
#endif
								((!m_pcEncCfg->getDisableIntraPUsInInterSlices()) && ((rpcBestCU->getCbf(0, COMPONENT_Y) != 0) ||
									((rpcBestCU->getCbf(0, COMPONENT_Cb) != 0) && (numberValidComponents > COMPONENT_Cb)) ||
									((rpcBestCU->getCbf(0, COMPONENT_Cr) != 0) && (numberValidComponents > COMPONENT_Cr)) // avoid very complex intra if it is unlikely
									)))
							{
								xCheckRDCostIntra(rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
								rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
								//if (uiDepth == sps.getLog2DiffMaxMinCodingBlockSize())
								//{
								//    if (rpcTempCU->getWidth(0) > (1 << sps.getQuadtreeTULog2MinSize()))
								//    {
								//        Double tmpIntraCost;
								//        xCheckRDCostIntra(rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug));
								//        intraCost = std::min(intraCost, tmpIntraCost);
								//        rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
								//    }
								//}
							}
						}

						double cost = rpcBestCU->getTotalCost();
						if (cost > 0)   // 确保cost有效
						{
							bool is_Intra;
							if (rpcBestCU->getPredictionMode(0) == 0) // ILR
								is_Intra = false;
							else if (rpcBestCU->getPredictionMode(0) == 1) // Intra
								is_Intra = true;

							if (uiDepth == 0)
							{
								curCosts0[px][py].first = (int)is_Intra;
								curCosts0[px][py].second = cost;
								cost0.second = cost;
							}
							else if (uiDepth == 1)
							{
								curCosts1[px1][py1].first = (int)is_Intra;
								curCosts1[px1][py1].second = cost;
								cost1[index[uiDepth]].second = cost;
							}
							else if (uiDepth == 2)
							{
								curCosts2[px2][py2].first = (int)is_Intra;
								curCosts2[px2][py2].second = cost;
								cost2[index[uiDepth]].second = cost;
							}
							else if (uiDepth == 3)
							{
								curCosts3[px3][py3].first = (int)is_Intra;
								curCosts3[px3][py3].second = cost;
							}
						}

						// 为啥先算pf[2]呢？本来是按顺序求的，但后来发现调整一下顺序速度会更快，因为相关性程度通常值较小，难以满足cond
						// 假如ITEW迭代结果和输入概率的顺序没关系，那么把序号改成顺序的应该对结果没啥影响，这里就不求证了
						if ((uiDepth == 1 || uiDepth == 2) && fr > 0 && px > 0 && px < n0 - 1 && py > 0 && py < m0 - 1)
						{
							pair<double, double> f3(-1, -1);
							if (curG[px - 1][py].isValid && curG[px][py - 1].isValid)
							{
								if (uiDepth == 1 && curmsg[px - 1][py].c1[2 + index[1] % 2].first == 0)
								{
									dg ug1 = curG[px - 1][py].lv1[1][index[1] % 2];
									f3.first = getCorrelation(ug1, cc);
								}
								else if (uiDepth == 2 && curmsg[px - 1][py].c2[12 + index[2] % 4].first == 0)
								{
									dg ug2 = curG[px - 1][py].lv2[3][index[2] % 4];
									f3.first = getCorrelation(ug2, cc);
								}
								if (f3.first > cond)
								{
									if (uiDepth == 1 && curmsg[px][py - 1].c1[1 + 2 * (int)(index[1] / 2)].first == 0)
									{
										dg lg1 = curG[px][py - 1].lv1[index[1] / 2][1];
										f3.second = getCorrelation(lg1, cc);
									}
									else if (uiDepth == 2 && curmsg[px][py - 1].c2[3 + 4 * (int)(index[2] / 4)].first == 0)
									{
										dg lg2 = curG[px][py - 1].lv2[index[2] / 4][3];
										f3.second = getCorrelation(lg2, cc);
									}
								}
								if (f3.second > cond)
									pf[2] = max(f3.first, f3.second);
								else
									pf[2] = -1;
							}
							else
								pf[2] = -1;
							if (pf[2] > cond)
							{
								pf[1] = hytest(rpcBestCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvBest[uiDepth], uiDepth);
								
								if (pf[1] > cond)
								{
									pf[0] = 0.499;
									if (uiDepth == 1)
									{
										vector<double> rec0;	// 终止
										vector<double> rec1;	// 分割
										int x, y;
										double e0 = 0, e1 = 0;
										if (premsg[px][py].isValid)
										{
											for (auto e : premsg[px][py].c1)
											{
												if (e.first == 0)
												{
													rec0.push_back(e.second);
													e0 += rec0.back();
												}
												else
												{
													rec1.push_back(e.second);
													e1 += rec1.back();
												}
											}
										}
										for (auto p : curPos)
										{
											x = px + p.first;
											y = py + p.second;
											if (curmsg[x][y].isValid)
											{
												for (auto e : curmsg[x][y].c1)
												{
													if (e.first == 0)
													{
														rec0.push_back(e.second);
														e0 += rec0.back();
													}
													else
													{
														rec1.push_back(e.second);
														e1 += rec1.back();
													}
												}
											}
										}
										if (rec0.size() < 1)
											pf[0] = 0;
										else if (rec1.size() < 1)
											pf[0] = 1;
										if (rec0.size() > 1 && rec1.size() > 1)
										{
											GMM gmm(rec0, rec1, e0, e1, cost);
											gmm.init();
											gmm.train();
											pf[0] = gmm.getResult();
										}
									}
									if (uiDepth == 2)
									{
										vector<double> rec0;	// 终止
										vector<double> rec1;	// 分割
										int x, y;
										double e0 = 0, e1 = 0;
										for (auto p : curPos2)
										{
											x = px + p.first;
											y = py + p.second;
											if (curmsg[x][y].isValid)
											{
												for (auto e : curmsg[x][y].c2)
												{
													if (e.first == 0)
													{
														rec0.push_back(e.second);
														e0 += rec0.back();
													}
													else
													{
														rec1.push_back(e.second);
														e1 += rec1.back();
													}
												}
											}
										}
										if (rec0.size() < 1)
											pf[0] = 0;
										else if (rec1.size() < 1)
											pf[0] = 1;
										if (rec0.size() > 1 && rec1.size() > 1)
										{
											GMM gmm(rec0, rec1, e0, e1, cost);
											gmm.init();
											gmm.train();
											pf[0] = gmm.getResult();
										}
									}
									if (pf[0] > cond)
									{

										ITEW ite(pf);
										vector<double> w = ite.getW();
										const double finalP = ite.getP();
										cease = (finalP > fcond) ? true : false;

									}
								}
							}
							
						}
						index[uiDepth]++;
					}
#endif
					// 基本层不管
					else
					{
						Double intraCost = 0.0;

						if ((rpcBestCU->getSlice()->getSliceType() == I_SLICE) ||
#if ENCODER_FAST_MODE
							rpcBestCU->getPredictionMode(0) == NUMBER_OF_PREDICTION_MODES || // if there is no valid inter prediction
							!testInter ||
#endif
							((!m_pcEncCfg->getDisableIntraPUsInInterSlices()) && ((rpcBestCU->getCbf(0, COMPONENT_Y) != 0) ||
								((rpcBestCU->getCbf(0, COMPONENT_Cb) != 0) && (numberValidComponents > COMPONENT_Cb)) ||
								((rpcBestCU->getCbf(0, COMPONENT_Cr) != 0) && (numberValidComponents > COMPONENT_Cr)) // avoid very complex intra if it is unlikely
								)))
						{
							xCheckRDCostIntra(rpcBestCU, rpcTempCU, intraCost, SIZE_2Nx2N DEBUG_STRING_PASS_INTO(sDebug));
							rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
							if (uiDepth == sps.getLog2DiffMaxMinCodingBlockSize())
							{
								if (rpcTempCU->getWidth(0) > (1 << sps.getQuadtreeTULog2MinSize()))
								{
									Double tmpIntraCost;
									xCheckRDCostIntra(rpcBestCU, rpcTempCU, tmpIntraCost, SIZE_NxN DEBUG_STRING_PASS_INTO(sDebug));
									intraCost = std::min(intraCost, tmpIntraCost);
									rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);
								}
							}
						}
					}
					if (bIsLosslessMode) // Restore loop variable if lossless mode was searched.
					{
						iQP = iMinQP;
					}
				}
			}

			if (rpcBestCU->getTotalCost() != MAX_DOUBLE)
			{
				m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
				m_pcEntropyCoder->resetBits();
				m_pcEntropyCoder->encodeSplitFlag(rpcBestCU, 0, uiDepth, true);
				rpcBestCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
				rpcBestCU->getTotalBins() += ((TEncBinCABAC*)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
				rpcBestCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcBestCU->getTotalBits(), rpcBestCU->getTotalDistortion());
				m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
			}
#if HIGHER_LAYER_IRAP_SKIP_FLAG
		}
#endif
	}

	// copy original YUV samples to PCM buffer
	if (rpcBestCU->getTotalCost() != MAX_DOUBLE && rpcBestCU->isLosslessCoded(0) && (rpcBestCU->getIPCMFlag(0) == false))
	{
		xFillPCMBuffer(rpcBestCU, m_ppcOrigYuv[uiDepth]);
	}

	if (uiDepth == pps.getMaxCuDQPDepth())
	{
		Int idQP = m_pcEncCfg->getMaxDeltaQP();
		iMinQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP - idQP);
		iMaxQP = Clip3(-sps.getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQP + idQP);
	}
	else if (uiDepth < pps.getMaxCuDQPDepth())
	{
		iMinQP = iBaseQP;
		iMaxQP = iBaseQP;
	}
	else
	{
		const Int iStartQP = rpcTempCU->getQP(0);
		iMinQP = iStartQP;
		iMaxQP = iStartQP;
	}

	if (m_pcEncCfg->getUseRateCtrl())
	{
		iMinQP = m_pcRateCtrl->getRCQP();
		iMaxQP = m_pcRateCtrl->getRCQP();
	}

	if (m_pcEncCfg->getCUTransquantBypassFlagForceValue())
	{
		iMaxQP = iMinQP; // If all TUs are forced into using transquant bypass, do not loop here.
	}

	const Bool bSubBranch = bBoundary || !(m_pcEncCfg->getUseEarlyCU() && rpcBestCU->getTotalCost() != MAX_DOUBLE && rpcBestCU->isSkipped(0));

	if (!cease && bSubBranch && uiDepth < sps.getLog2DiffMaxMinCodingBlockSize() && (!getFastDeltaQp() || uiWidth > fastDeltaQPCuMaxSize || bBoundary))
	{
		// further split
		for (Int iQP = iMinQP; iQP <= iMaxQP; iQP++)
		{
			const Bool bIsLosslessMode = false; // False at this level. Next level down may set it to true.

			rpcTempCU->initEstData(uiDepth, iQP, bIsLosslessMode);

			UChar uhNextDepth = uiDepth + 1;
			TComDataCU* pcSubBestPartCU = m_ppcBestCU[uhNextDepth];
			TComDataCU* pcSubTempPartCU = m_ppcTempCU[uhNextDepth];
			DEBUG_STRING_NEW(sTempDebug)

				for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++)
				{
					pcSubBestPartCU->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP); // clear sub partition datas or init.
					pcSubTempPartCU->initSubCU(rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP); // clear sub partition datas or init.

					if ((pcSubBestPartCU->getCUPelX() < sps.getPicWidthInLumaSamples()) && (pcSubBestPartCU->getCUPelY() < sps.getPicHeightInLumaSamples()))
					{
						if (0 == uiPartUnitIdx) //initialize RD with previous depth buffer
						{
							m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
						}
						else
						{
							m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]);
						}

#if AMP_ENC_SPEEDUP
						DEBUG_STRING_NEW(sChild)
							if (!(rpcBestCU->getTotalCost() != MAX_DOUBLE && rpcBestCU->isInter(0)))
							{
								xCompressCU(pcSubBestPartCU, pcSubTempPartCU, uhNextDepth DEBUG_STRING_PASS_INTO(sChild), NUMBER_OF_PART_SIZES);
							}
							else
							{

								xCompressCU(pcSubBestPartCU, pcSubTempPartCU, uhNextDepth DEBUG_STRING_PASS_INTO(sChild), rpcBestCU->getPartitionSize(0));
							}
						DEBUG_STRING_APPEND(sTempDebug, sChild)
#else
						xCompressCU(pcSubBestPartCU, pcSubTempPartCU, uhNextDepth);
#endif

						rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth); // Keep best part data to current temporary data.
						xCopyYuv2Tmp(pcSubBestPartCU->getTotalNumPart() * uiPartUnitIdx, uhNextDepth);
					}
					else
					{
						pcSubBestPartCU->copyToPic(uhNextDepth);
						rpcTempCU->copyPartFrom(pcSubBestPartCU, uiPartUnitIdx, uhNextDepth);
					}
				}

			m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]);
			if (!bBoundary)
			{
				m_pcEntropyCoder->resetBits();
				m_pcEntropyCoder->encodeSplitFlag(rpcTempCU, 0, uiDepth, true);

				rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
				rpcTempCU->getTotalBins() += ((TEncBinCABAC*)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
			}
			rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());

			if (uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
			{
				Bool hasResidual = false;
				for (UInt uiBlkIdx = 0; uiBlkIdx < rpcTempCU->getTotalNumPart(); uiBlkIdx++)
				{
					if ((rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Y) || (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cb) && (numberValidComponents > COMPONENT_Cb)) || (rpcTempCU->getCbf(uiBlkIdx, COMPONENT_Cr) && (numberValidComponents > COMPONENT_Cr))))
					{
						hasResidual = true;
						break;
					}
				}

				if (hasResidual)
				{
					m_pcEntropyCoder->resetBits();
					m_pcEntropyCoder->encodeQP(rpcTempCU, 0, false);
					rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
					rpcTempCU->getTotalBins() += ((TEncBinCABAC*)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
					rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());

					Bool foundNonZeroCbf = false;
					rpcTempCU->setQPSubCUs(rpcTempCU->getRefQP(0), 0, uiDepth, foundNonZeroCbf);
					assert(foundNonZeroCbf);
				}
				else
				{
					rpcTempCU->setQPSubParts(rpcTempCU->getRefQP(0), 0, uiDepth); // set QP to default QP
				}
			}

			m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

			// If the configuration being tested exceeds the maximum number of bytes for a slice / slice-segment, then
			// a proper RD evaluation cannot be performed. Therefore, termination of the
			// slice/slice-segment must be made prior to this CTU.
			// This can be achieved by forcing the decision to be that of the rpcTempCU.
			// The exception is each slice / slice-segment must have at least one CTU.
			if (rpcBestCU->getTotalCost() != MAX_DOUBLE)
			{
				const Bool isEndOfSlice = pcSlice->getSliceMode() == FIXED_NUMBER_OF_BYTES && ((pcSlice->getSliceBits() + rpcBestCU->getTotalBits()) > pcSlice->getSliceArgument() << 3) && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceCurStartCtuTsAddr()) && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
				const Bool isEndOfSliceSegment = pcSlice->getSliceSegmentMode() == FIXED_NUMBER_OF_BYTES && ((pcSlice->getSliceSegmentBits() + rpcBestCU->getTotalBits()) > pcSlice->getSliceSegmentArgument() << 3) && rpcBestCU->getCtuRsAddr() != pcPic->getPicSym()->getCtuTsToRsAddrMap(pcSlice->getSliceSegmentCurStartCtuTsAddr());
				// Do not need to check slice condition for slice-segment since a slice-segment is a subset of a slice.
				if (isEndOfSlice || isEndOfSliceSegment)
				{
					rpcBestCU->getTotalCost() = MAX_DOUBLE;
				}
			}

			xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTempDebug) DEBUG_STRING_PASS_INTO(false)); // RD compare current larger prediction
																																						   // with sub partitioned prediction.
		}
	}

	DEBUG_STRING_APPEND(sDebug_, sDebug);

	rpcBestCU->copyToPic(uiDepth); // Copy Best data to Picture for next partition prediction.

	xCopyYuv2Pic(rpcBestCU->getPic(), rpcBestCU->getCtuRsAddr(), rpcBestCU->getZorderIdxInCtu(), uiDepth, uiDepth); // Copy Yuv data to picture Yuv
	if (bBoundary)
	{
		return;
	}

	// Assert if Best prediction mode is NONE
	// Selected mode's RD-cost must be not MAX_DOUBLE.
	assert(rpcBestCU->getPartitionSize(0) != NUMBER_OF_PART_SIZES);
	assert(rpcBestCU->getPredictionMode(0) != NUMBER_OF_PREDICTION_MODES);
	assert(rpcBestCU->getTotalCost() != MAX_DOUBLE);
}

/** finish encoding a cu and handle end-of-slice conditions
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TEncCu::finishCU(TComDataCU* pcCU, UInt uiAbsPartIdx)
{
	TComPic* pcPic = pcCU->getPic();
	TComSlice* pcSlice = pcCU->getPic()->getSlice(pcCU->getPic()->getCurrSliceIdx());

	//Calculate end address
	const Int currentCTUTsAddr = pcPic->getPicSym()->getCtuRsToTsAddrMap(pcCU->getCtuRsAddr());
	const Bool isLastSubCUOfCtu = pcCU->isLastSubCUOfCtu(uiAbsPartIdx);
	if (isLastSubCUOfCtu)
	{
		// The 1-terminating bit is added to all streams, so don't add it here when it's 1.
		// i.e. when the slice segment CurEnd CTU address is the current CTU address+1.
		if (pcSlice->getSliceSegmentCurEndCtuTsAddr() != currentCTUTsAddr + 1)
		{
			m_pcEntropyCoder->encodeTerminatingBit(0);
		}
	}
}

/** Compute QP for each CU
 * \param pcCU Target CU
 * \param uiDepth CU depth
 * \returns quantization parameter
 */
Int TEncCu::xComputeQP(TComDataCU* pcCU, UInt uiDepth)
{
	Int iBaseQp = pcCU->getSlice()->getSliceQp();
	Int iQpOffset = 0;
	if (m_pcEncCfg->getUseAdaptiveQP())
	{
		TEncPic* pcEPic = dynamic_cast<TEncPic*>(pcCU->getPic());
		UInt uiAQDepth = min(uiDepth, pcEPic->getMaxAQDepth() - 1);
		TEncPicQPAdaptationLayer* pcAQLayer = pcEPic->getAQLayer(uiAQDepth);
		UInt uiAQUPosX = pcCU->getCUPelX() / pcAQLayer->getAQPartWidth();
		UInt uiAQUPosY = pcCU->getCUPelY() / pcAQLayer->getAQPartHeight();
		UInt uiAQUStride = pcAQLayer->getAQPartStride();
		TEncQPAdaptationUnit* acAQU = pcAQLayer->getQPAdaptationUnit();

		Double dMaxQScale = pow(2.0, m_pcEncCfg->getQPAdaptationRange() / 6.0);
		Double dAvgAct = pcAQLayer->getAvgActivity();
		Double dCUAct = acAQU[uiAQUPosY * uiAQUStride + uiAQUPosX].getActivity();
		Double dNormAct = (dMaxQScale * dCUAct + dAvgAct) / (dCUAct + dMaxQScale * dAvgAct);
		Double dQpOffset = log(dNormAct) / log(2.0) * 6.0;
		iQpOffset = Int(floor(dQpOffset + 0.49999));
	}
	return Clip3(-pcCU->getSlice()->getSPS()->getQpBDOffset(CHANNEL_TYPE_LUMA), MAX_QP, iBaseQp + iQpOffset);
}

/** encode a CU block recursively
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth
 * \returns Void
 */
Void TEncCu::xEncodeCU(TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth)
{
	TComPic* const pcPic = pcCU->getPic();
	TComSlice* const pcSlice = pcCU->getSlice();
	const TComSPS& sps = *(pcSlice->getSPS());
	const TComPPS& pps = *(pcSlice->getPPS());

	const UInt maxCUWidth = sps.getMaxCUWidth();
	const UInt maxCUHeight = sps.getMaxCUHeight();

	Bool bBoundary = false;
	UInt uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
	const UInt uiRPelX = uiLPelX + (maxCUWidth >> uiDepth) - 1;
	UInt uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];
	const UInt uiBPelY = uiTPelY + (maxCUHeight >> uiDepth) - 1;

#if HIGHER_LAYER_IRAP_SKIP_FLAG
	if (m_pcEncCfg->getSkipPictureAtArcSwitch() && m_pcEncCfg->getAdaptiveResolutionChange() > 0 && pcSlice->getLayerId() == 1 && pcSlice->getPOC() == m_pcEncCfg->getAdaptiveResolutionChange())
	{
		pcCU->setSkipFlagSubParts(true, uiAbsPartIdx, uiDepth);
	}
#endif

	if ((uiRPelX < sps.getPicWidthInLumaSamples()) && (uiBPelY < sps.getPicHeightInLumaSamples()))
	{
		m_pcEntropyCoder->encodeSplitFlag(pcCU, uiAbsPartIdx, uiDepth);
	}
	else
	{
		bBoundary = true;
	}

	if (((uiDepth < pcCU->getDepth(uiAbsPartIdx)) && (uiDepth < sps.getLog2DiffMaxMinCodingBlockSize())) || bBoundary)
	{
		UInt uiQNumParts = (pcPic->getNumPartitionsInCtu() >> (uiDepth << 1)) >> 2;
		if (uiDepth == pps.getMaxCuDQPDepth() && pps.getUseDQP())
		{
			setdQPFlag(true);
		}

		if (uiDepth == pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcSlice->getUseChromaQpAdj())
		{
			setCodeChromaQpAdjFlag(true);
		}

		for (UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx += uiQNumParts)
		{
			uiLPelX = pcCU->getCUPelX() + g_auiRasterToPelX[g_auiZscanToRaster[uiAbsPartIdx]];
			uiTPelY = pcCU->getCUPelY() + g_auiRasterToPelY[g_auiZscanToRaster[uiAbsPartIdx]];

			if ((uiLPelX < sps.getPicWidthInLumaSamples()) && (uiTPelY < sps.getPicHeightInLumaSamples()))
			{
				xEncodeCU(pcCU, uiAbsPartIdx, uiDepth + 1);
			}
		}
		return;
	}

	if (uiDepth <= pps.getMaxCuDQPDepth() && pps.getUseDQP())
	{
		setdQPFlag(true);
	}

	if (uiDepth <= pps.getPpsRangeExtension().getDiffCuChromaQpOffsetDepth() && pcSlice->getUseChromaQpAdj())
	{
		setCodeChromaQpAdjFlag(true);
	}

	if (pps.getTransquantBypassEnableFlag())
	{
		m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, uiAbsPartIdx);
	}

	if (!pcSlice->isIntra())
	{
		m_pcEntropyCoder->encodeSkipFlag(pcCU, uiAbsPartIdx);
	}

	if (pcCU->isSkipped(uiAbsPartIdx))
	{
		m_pcEntropyCoder->encodeMergeIndex(pcCU, uiAbsPartIdx);
		finishCU(pcCU, uiAbsPartIdx);
		return;
	}

	m_pcEntropyCoder->encodePredMode(pcCU, uiAbsPartIdx);
	m_pcEntropyCoder->encodePartSize(pcCU, uiAbsPartIdx, uiDepth);

	if (pcCU->isIntra(uiAbsPartIdx) && pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N)
	{
		m_pcEntropyCoder->encodeIPCMInfo(pcCU, uiAbsPartIdx);

		if (pcCU->getIPCMFlag(uiAbsPartIdx))
		{
			// Encode slice finish
			finishCU(pcCU, uiAbsPartIdx);
			return;
		}
	}

	// prediction Info ( Intra : direction mode, Inter : Mv, reference idx )
	m_pcEntropyCoder->encodePredInfo(pcCU, uiAbsPartIdx);

	// Encode Coefficients
	Bool bCodeDQP = getdQPFlag();
	Bool codeChromaQpAdj = getCodeChromaQpAdjFlag();
	m_pcEntropyCoder->encodeCoeff(pcCU, uiAbsPartIdx, uiDepth, bCodeDQP, codeChromaQpAdj);
	setCodeChromaQpAdjFlag(codeChromaQpAdj);
	setdQPFlag(bCodeDQP);

	// --- write terminating bit ---
	finishCU(pcCU, uiAbsPartIdx);
}

Int xCalcHADs8x8_ISlice(Pel* piOrg, Int iStrideOrg)
{
	Int k, i, j, jj;
	Int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

	for (k = 0; k < 64; k += 8)
	{
		diff[k + 0] = piOrg[0];
		diff[k + 1] = piOrg[1];
		diff[k + 2] = piOrg[2];
		diff[k + 3] = piOrg[3];
		diff[k + 4] = piOrg[4];
		diff[k + 5] = piOrg[5];
		diff[k + 6] = piOrg[6];
		diff[k + 7] = piOrg[7];

		piOrg += iStrideOrg;
	}

	//horizontal
	for (j = 0; j < 8; j++)
	{
		jj = j << 3;
		m2[j][0] = diff[jj] + diff[jj + 4];
		m2[j][1] = diff[jj + 1] + diff[jj + 5];
		m2[j][2] = diff[jj + 2] + diff[jj + 6];
		m2[j][3] = diff[jj + 3] + diff[jj + 7];
		m2[j][4] = diff[jj] - diff[jj + 4];
		m2[j][5] = diff[jj + 1] - diff[jj + 5];
		m2[j][6] = diff[jj + 2] - diff[jj + 6];
		m2[j][7] = diff[jj + 3] - diff[jj + 7];

		m1[j][0] = m2[j][0] + m2[j][2];
		m1[j][1] = m2[j][1] + m2[j][3];
		m1[j][2] = m2[j][0] - m2[j][2];
		m1[j][3] = m2[j][1] - m2[j][3];
		m1[j][4] = m2[j][4] + m2[j][6];
		m1[j][5] = m2[j][5] + m2[j][7];
		m1[j][6] = m2[j][4] - m2[j][6];
		m1[j][7] = m2[j][5] - m2[j][7];

		m2[j][0] = m1[j][0] + m1[j][1];
		m2[j][1] = m1[j][0] - m1[j][1];
		m2[j][2] = m1[j][2] + m1[j][3];
		m2[j][3] = m1[j][2] - m1[j][3];
		m2[j][4] = m1[j][4] + m1[j][5];
		m2[j][5] = m1[j][4] - m1[j][5];
		m2[j][6] = m1[j][6] + m1[j][7];
		m2[j][7] = m1[j][6] - m1[j][7];
	}

	//vertical
	for (i = 0; i < 8; i++)
	{
		m3[0][i] = m2[0][i] + m2[4][i];
		m3[1][i] = m2[1][i] + m2[5][i];
		m3[2][i] = m2[2][i] + m2[6][i];
		m3[3][i] = m2[3][i] + m2[7][i];
		m3[4][i] = m2[0][i] - m2[4][i];
		m3[5][i] = m2[1][i] - m2[5][i];
		m3[6][i] = m2[2][i] - m2[6][i];
		m3[7][i] = m2[3][i] - m2[7][i];

		m1[0][i] = m3[0][i] + m3[2][i];
		m1[1][i] = m3[1][i] + m3[3][i];
		m1[2][i] = m3[0][i] - m3[2][i];
		m1[3][i] = m3[1][i] - m3[3][i];
		m1[4][i] = m3[4][i] + m3[6][i];
		m1[5][i] = m3[5][i] + m3[7][i];
		m1[6][i] = m3[4][i] - m3[6][i];
		m1[7][i] = m3[5][i] - m3[7][i];

		m2[0][i] = m1[0][i] + m1[1][i];
		m2[1][i] = m1[0][i] - m1[1][i];
		m2[2][i] = m1[2][i] + m1[3][i];
		m2[3][i] = m1[2][i] - m1[3][i];
		m2[4][i] = m1[4][i] + m1[5][i];
		m2[5][i] = m1[4][i] - m1[5][i];
		m2[6][i] = m1[6][i] + m1[7][i];
		m2[7][i] = m1[6][i] - m1[7][i];
	}

	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			iSumHad += abs(m2[i][j]);
		}
	}
	iSumHad -= abs(m2[0][0]);
	iSumHad = (iSumHad + 2) >> 2;
	return (iSumHad);
}

Int TEncCu::updateCtuDataISlice(TComDataCU* pCtu, Int width, Int height)
{
	Int xBl, yBl;
	const Int iBlkSize = 8;

	Pel* pOrgInit = pCtu->getPic()->getPicYuvOrg()->getAddr(COMPONENT_Y, pCtu->getCtuRsAddr(), 0);
	Int iStrideOrig = pCtu->getPic()->getPicYuvOrg()->getStride(COMPONENT_Y);
	Pel* pOrg;

	Int iSumHad = 0;
	for (yBl = 0; (yBl + iBlkSize) <= height; yBl += iBlkSize)
	{
		for (xBl = 0; (xBl + iBlkSize) <= width; xBl += iBlkSize)
		{
			pOrg = pOrgInit + iStrideOrig * yBl + xBl;
			iSumHad += xCalcHADs8x8_ISlice(pOrg, iStrideOrig);
		}
	}
	return (iSumHad);
}

/** check RD costs for a CU block encoded with merge
 * \param rpcBestCU
 * \param rpcTempCU
 * \param earlyDetectionSkipMode
 */
#if HIGHER_LAYER_IRAP_SKIP_FLAG
Void TEncCu::xCheckRDCostMerge2Nx2N(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU DEBUG_STRING_FN_DECLARE(sDebug), Bool* earlyDetectionSkipMode, Bool bUseSkip)
#else
Void TEncCu::xCheckRDCostMerge2Nx2N(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU DEBUG_STRING_FN_DECLARE(sDebug), Bool* earlyDetectionSkipMode)
#endif
{
	assert(rpcTempCU->getSlice()->getSliceType() != I_SLICE);
	if (getFastDeltaQp())
	{
		return; // never check merge in fast deltaqp mode
	}
	TComMvField cMvFieldNeighbours[2 * MRG_MAX_NUM_CANDS]; // double length for mv of both lists
	UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
	Int numValidMergeCand = 0;
	const Bool bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);

	for (UInt ui = 0; ui < rpcTempCU->getSlice()->getMaxNumMergeCand(); ++ui)
	{
		uhInterDirNeighbours[ui] = 0;
	}
	UChar uhDepth = rpcTempCU->getDepth(0);
	rpcTempCU->setPartSizeSubParts(SIZE_2Nx2N, 0, uhDepth); // interprets depth relative to CTU level
	rpcTempCU->getInterMergeCandidates(0, 0, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand);

	Int mergeCandBuffer[MRG_MAX_NUM_CANDS];
	for (UInt ui = 0; ui < numValidMergeCand; ++ui)
	{
		mergeCandBuffer[ui] = 0;
	}

	Bool bestIsSkip = false;

	UInt iteration;
	if (rpcTempCU->isLosslessCoded(0))
	{
		iteration = 1;
	}
	else
	{
		iteration = 2;
	}
	DEBUG_STRING_NEW(bestStr)

#if HIGHER_LAYER_IRAP_SKIP_FLAG
		for (UInt uiNoResidual = bUseSkip ? 1 : 0; uiNoResidual < iteration; ++uiNoResidual)
#else
		for (UInt uiNoResidual = 0; uiNoResidual < iteration; ++uiNoResidual)
#endif
		{
			for (UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand)
			{
#if REF_IDX_ME_ZEROMV
				Bool bZeroMVILR = rpcTempCU->xCheckZeroMVILRMerge(uhInterDirNeighbours[uiMergeCand], cMvFieldNeighbours[0 + 2 * uiMergeCand], cMvFieldNeighbours[1 + 2 * uiMergeCand]);
				if (bZeroMVILR)
				{
#endif
#if N0383_IL_CONSTRAINED_TILE_SETS_SEI
					if (!(rpcTempCU->isInterLayerReference(uhInterDirNeighbours[uiMergeCand], cMvFieldNeighbours[0 + 2 * uiMergeCand], cMvFieldNeighbours[1 + 2 * uiMergeCand]) && m_disableILP))
					{
#endif
						if (!(uiNoResidual == 1 && mergeCandBuffer[uiMergeCand] == 1))
						{
							if (!(bestIsSkip && uiNoResidual == 0))
							{
								DEBUG_STRING_NEW(tmpStr)
									// set MC parameters
									rpcTempCU->setPredModeSubParts(MODE_INTER, 0, uhDepth); // interprets depth relative to CTU level
								rpcTempCU->setCUTransquantBypassSubParts(bTransquantBypassFlag, 0, uhDepth);
								rpcTempCU->setChromaQpAdjSubParts(bTransquantBypassFlag ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth);
								rpcTempCU->setPartSizeSubParts(SIZE_2Nx2N, 0, uhDepth);                                                            // interprets depth relative to CTU level
								rpcTempCU->setMergeFlagSubParts(true, 0, 0, uhDepth);                                                              // interprets depth relative to CTU level
								rpcTempCU->setMergeIndexSubParts(uiMergeCand, 0, 0, uhDepth);                                                      // interprets depth relative to CTU level
								rpcTempCU->setInterDirSubParts(uhInterDirNeighbours[uiMergeCand], 0, 0, uhDepth);                                  // interprets depth relative to CTU level
								rpcTempCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField(cMvFieldNeighbours[0 + 2 * uiMergeCand], SIZE_2Nx2N, 0, 0); // interprets depth relative to rpcTempCU level
								rpcTempCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField(cMvFieldNeighbours[1 + 2 * uiMergeCand], SIZE_2Nx2N, 0, 0); // interprets depth relative to rpcTempCU level

								// do MC
								m_pcPredSearch->motionCompensation(rpcTempCU, m_ppcPredYuvTemp[uhDepth]);
								// estimate residual and encode everything
								m_pcPredSearch->encodeResAndCalcRdInterCU(rpcTempCU,
									m_ppcOrigYuv[uhDepth],
									m_ppcPredYuvTemp[uhDepth],
									m_ppcResiYuvTemp[uhDepth],
									m_ppcResiYuvBest[uhDepth],
									m_ppcRecoYuvTemp[uhDepth],
									(uiNoResidual != 0) DEBUG_STRING_PASS_INTO(tmpStr));

#if DEBUG_STRING
								DebugInterPredResiReco(tmpStr, *(m_ppcPredYuvTemp[uhDepth]), *(m_ppcResiYuvBest[uhDepth]), *(m_ppcRecoYuvTemp[uhDepth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif

								if ((uiNoResidual == 0) && (rpcTempCU->getQtRootCbf(0) == 0))
								{
									// If no residual when allowing for one, then set mark to not try case where residual is forced to 0
									mergeCandBuffer[uiMergeCand] = 1;
								}

								Int orgQP = rpcTempCU->getQP(0);
								xCheckDQP(rpcTempCU);
								xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(bestStr) DEBUG_STRING_PASS_INTO(tmpStr));

								rpcTempCU->initEstData(uhDepth, orgQP, bTransquantBypassFlag);

								if (m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip)
								{
									bestIsSkip = rpcBestCU->getQtRootCbf(0) == 0;
								}
							}
						}
#if N0383_IL_CONSTRAINED_TILE_SETS_SEI
					}
#endif
#if REF_IDX_ME_ZEROMV
				}
#endif
			}

			if (uiNoResidual == 0 && m_pcEncCfg->getUseEarlySkipDetection())
			{
				if (rpcBestCU->getQtRootCbf(0) == 0)
				{
					if (rpcBestCU->getMergeFlag(0))
					{
						*earlyDetectionSkipMode = true;
					}
					else if (m_pcEncCfg->getMotionEstimationSearchMethod() != MESEARCH_SELECTIVE)
					{
						Int absoulte_MV = 0;
						for (UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
						{
							if (rpcBestCU->getSlice()->getNumRefIdx(RefPicList(uiRefListIdx)) > 0)
							{
								TComCUMvField* pcCUMvField = rpcBestCU->getCUMvField(RefPicList(uiRefListIdx));
								Int iHor = pcCUMvField->getMvd(0).getAbsHor();
								Int iVer = pcCUMvField->getMvd(0).getAbsVer();
								absoulte_MV += iHor + iVer;
							}
						}

						if (absoulte_MV == 0)
						{
							*earlyDetectionSkipMode = true;
						}
					}
				}
			}
		}
	DEBUG_STRING_APPEND(sDebug, bestStr)
}

#if AMP_MRG
Void TEncCu::xCheckRDCostInter(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize DEBUG_STRING_FN_DECLARE(sDebug), Bool bUseMRG)
#else
Void TEncCu::xCheckRDCostInter(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize)
#endif
{
	DEBUG_STRING_NEW(sTest)

		if (getFastDeltaQp())
		{
			const TComSPS& sps = *(rpcTempCU->getSlice()->getSPS());
			const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight() >> (sps.getLog2DiffMaxMinCodingBlockSize()), sps.getMaxCUHeight(), 32u);
			if (ePartSize != SIZE_2Nx2N || rpcTempCU->getWidth(0) > fastDeltaQPCuMaxSize)
			{
				return; // only check necessary 2Nx2N Inter in fast deltaqp mode
			}
		}

	// prior to this, rpcTempCU will have just been reset using rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
	UChar uhDepth = rpcTempCU->getDepth(0);

	rpcTempCU->setPartSizeSubParts(ePartSize, 0, uhDepth);
	rpcTempCU->setPredModeSubParts(MODE_INTER, 0, uhDepth);
	rpcTempCU->setChromaQpAdjSubParts(rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uhDepth);

#if SVC_EXTENSION
#if AMP_MRG
	rpcTempCU->setMergeAMP(true);
	Bool ret = m_pcPredSearch->predInterSearch(rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth] DEBUG_STRING_PASS_INTO(sTest), false, bUseMRG);
#else
	Bool ret = m_pcPredSearch->predInterSearch(rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth]);
#endif

	if (!ret)
	{
		return;
	}
#else
#if AMP_MRG
	rpcTempCU->setMergeAMP(true);
	m_pcPredSearch->predInterSearch(rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth] DEBUG_STRING_PASS_INTO(sTest), false, bUseMRG);
#else
	m_pcPredSearch->predInterSearch(rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth]);
#endif
#endif

#if AMP_MRG
	if (!rpcTempCU->getMergeAMP())
	{
		return;
	}
#endif

	m_pcPredSearch->encodeResAndCalcRdInterCU(rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false DEBUG_STRING_PASS_INTO(sTest));
	rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());

#if DEBUG_STRING
	DebugInterPredResiReco(sTest, *(m_ppcPredYuvTemp[uhDepth]), *(m_ppcResiYuvBest[uhDepth]), *(m_ppcRecoYuvTemp[uhDepth]), DebugStringGetPredModeMask(rpcTempCU->getPredictionMode(0)));
#endif

	xCheckDQP(rpcTempCU);
	xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
}

Void TEncCu::xCheckRDCostIntra(TComDataCU*& rpcBestCU,
	TComDataCU*& rpcTempCU,
	Double& cost,
	PartSize eSize
	DEBUG_STRING_FN_DECLARE(sDebug))
{
	DEBUG_STRING_NEW(sTest)

		if (getFastDeltaQp())
		{
			const TComSPS& sps = *(rpcTempCU->getSlice()->getSPS());
			const UInt fastDeltaQPCuMaxSize = Clip3(sps.getMaxCUHeight() >> (sps.getLog2DiffMaxMinCodingBlockSize()), sps.getMaxCUHeight(), 32u);
			if (rpcTempCU->getWidth(0) > fastDeltaQPCuMaxSize)
			{
				return; // only check necessary 2Nx2N Intra in fast deltaqp mode
			}
		}

	UInt uiDepth = rpcTempCU->getDepth(0);

	rpcTempCU->setSkipFlagSubParts(false, 0, uiDepth);

	rpcTempCU->setPartSizeSubParts(eSize, 0, uiDepth);
	rpcTempCU->setPredModeSubParts(MODE_INTRA, 0, uiDepth);
	rpcTempCU->setChromaQpAdjSubParts(rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uiDepth);

	Pel resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];

	m_pcPredSearch->estIntraPredLumaQT(rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma DEBUG_STRING_PASS_INTO(sTest));

	m_ppcRecoYuvTemp[uiDepth]->copyToPicComponent(COMPONENT_Y, rpcTempCU->getPic()->getPicYuvRec(), rpcTempCU->getCtuRsAddr(), rpcTempCU->getZorderIdxInCtu());

	if (rpcBestCU->getPic()->getChromaFormat() != CHROMA_400)
	{
		m_pcPredSearch->estIntraPredChromaQT(rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], resiLuma DEBUG_STRING_PASS_INTO(sTest));
	}

	m_pcEntropyCoder->resetBits();

	if (rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
	{
		m_pcEntropyCoder->encodeCUTransquantBypassFlag(rpcTempCU, 0, true);
	}

	m_pcEntropyCoder->encodeSkipFlag(rpcTempCU, 0, true);
	m_pcEntropyCoder->encodePredMode(rpcTempCU, 0, true);
	m_pcEntropyCoder->encodePartSize(rpcTempCU, 0, uiDepth, true);
	m_pcEntropyCoder->encodePredInfo(rpcTempCU, 0);
	m_pcEntropyCoder->encodeIPCMInfo(rpcTempCU, 0, true);

	// Encode Coefficients
	Bool bCodeDQP = getdQPFlag();
	Bool codeChromaQpAdjFlag = getCodeChromaQpAdjFlag();
	m_pcEntropyCoder->encodeCoeff(rpcTempCU, 0, uiDepth, bCodeDQP, codeChromaQpAdjFlag);
	setCodeChromaQpAdjFlag(codeChromaQpAdjFlag);
	setdQPFlag(bCodeDQP);

	m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

	rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
	rpcTempCU->getTotalBins() = ((TEncBinCABAC*)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
	rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());

	xCheckDQP(rpcTempCU);

	cost = rpcTempCU->getTotalCost();

	xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(sDebug) DEBUG_STRING_PASS_INTO(sTest));
}

/** Check R-D costs for a CU with PCM mode.
 * \param rpcBestCU pointer to best mode CU data structure
 * \param rpcTempCU pointer to testing mode CU data structure
 * \returns Void
 *
 * \note Current PCM implementation encodes sample values in a lossless way. The distortion of PCM mode CUs are zero. PCM mode is selected if the best mode yields bits greater than that of PCM mode.
 */
Void TEncCu::xCheckIntraPCM(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU)
{
	if (getFastDeltaQp())
	{
		const TComSPS& sps = *(rpcTempCU->getSlice()->getSPS());
		const UInt fastDeltaQPCuMaxPCMSize = Clip3((UInt)1 << sps.getPCMLog2MinSize(), (UInt)1 << sps.getPCMLog2MaxSize(), 32u);
		if (rpcTempCU->getWidth(0) > fastDeltaQPCuMaxPCMSize)
		{
			return; // only check necessary PCM in fast deltaqp mode
		}
	}

	UInt uiDepth = rpcTempCU->getDepth(0);

	rpcTempCU->setSkipFlagSubParts(false, 0, uiDepth);

	rpcTempCU->setIPCMFlag(0, true);
	rpcTempCU->setIPCMFlagSubParts(true, 0, rpcTempCU->getDepth(0));
	rpcTempCU->setPartSizeSubParts(SIZE_2Nx2N, 0, uiDepth);
	rpcTempCU->setPredModeSubParts(MODE_INTRA, 0, uiDepth);
	rpcTempCU->setTrIdxSubParts(0, 0, uiDepth);
	rpcTempCU->setChromaQpAdjSubParts(rpcTempCU->getCUTransquantBypass(0) ? 0 : m_cuChromaQpOffsetIdxPlus1, 0, uiDepth);

	m_pcPredSearch->IPCMSearch(rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth]);

	m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);

	m_pcEntropyCoder->resetBits();

	if (rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
	{
		m_pcEntropyCoder->encodeCUTransquantBypassFlag(rpcTempCU, 0, true);
	}

	m_pcEntropyCoder->encodeSkipFlag(rpcTempCU, 0, true);
	m_pcEntropyCoder->encodePredMode(rpcTempCU, 0, true);
	m_pcEntropyCoder->encodePartSize(rpcTempCU, 0, uiDepth, true);
	m_pcEntropyCoder->encodeIPCMInfo(rpcTempCU, 0, true);

	m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

	rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
	rpcTempCU->getTotalBins() = ((TEncBinCABAC*)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
	rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());

	xCheckDQP(rpcTempCU);
	DEBUG_STRING_NEW(a)
		DEBUG_STRING_NEW(b)
		xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth DEBUG_STRING_PASS_INTO(a) DEBUG_STRING_PASS_INTO(b));
}

/** check whether current try is the best with identifying the depth of current try
 * \param rpcBestCU
 * \param rpcTempCU
 * \param uiDepth
 */
Void TEncCu::xCheckBestMode(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth DEBUG_STRING_FN_DECLARE(sParent) DEBUG_STRING_FN_DECLARE(sTest) DEBUG_STRING_PASS_INTO(Bool bAddSizeInfo))
{
	if (rpcTempCU->getTotalCost() < rpcBestCU->getTotalCost())
	{
		TComYuv* pcYuv;
		// Change Information data
		TComDataCU* pcCU = rpcBestCU;
		rpcBestCU = rpcTempCU;
		rpcTempCU = pcCU;

		// Change Prediction data
		pcYuv = m_ppcPredYuvBest[uiDepth];
		m_ppcPredYuvBest[uiDepth] = m_ppcPredYuvTemp[uiDepth];
		m_ppcPredYuvTemp[uiDepth] = pcYuv;

		// Change Reconstruction data
		pcYuv = m_ppcRecoYuvBest[uiDepth];
		m_ppcRecoYuvBest[uiDepth] = m_ppcRecoYuvTemp[uiDepth];
		m_ppcRecoYuvTemp[uiDepth] = pcYuv;

		pcYuv = NULL;
		pcCU = NULL;

		// store temp best CI for next CU coding
		m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]->store(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);

#if DEBUG_STRING
		DEBUG_STRING_SWAP(sParent, sTest)
			const PredMode predMode = rpcBestCU->getPredictionMode(0);
		if ((DebugOptionList::DebugString_Structure.getInt() & DebugStringGetPredModeMask(predMode)) && bAddSizeInfo)
		{
			std::stringstream ss(stringstream::out);
			ss << "###: " << (predMode == MODE_INTRA ? "Intra   " : "Inter   ") << partSizeToString[rpcBestCU->getPartitionSize(0)] << " CU at " << rpcBestCU->getCUPelX() << ", " << rpcBestCU->getCUPelY() << " width=" << UInt(rpcBestCU->getWidth(0)) << std::endl;
			sParent += ss.str();
		}
#endif
	}
}

Void TEncCu::xCheckDQP(TComDataCU* pcCU)
{
	UInt uiDepth = pcCU->getDepth(0);

	const TComPPS& pps = *(pcCU->getSlice()->getPPS());
	if (pps.getUseDQP() && uiDepth <= pps.getMaxCuDQPDepth())
	{
		if (pcCU->getQtRootCbf(0))
		{
			m_pcEntropyCoder->resetBits();
			m_pcEntropyCoder->encodeQP(pcCU, 0, false);
			pcCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
			pcCU->getTotalBins() += ((TEncBinCABAC*)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
			pcCU->getTotalCost() = m_pcRdCost->calcRdCost(pcCU->getTotalBits(), pcCU->getTotalDistortion());
		}
		else
		{
			pcCU->setQPSubParts(pcCU->getRefQP(0), 0, uiDepth); // set QP to default QP
		}
	}
}

Void TEncCu::xCopyAMVPInfo(AMVPInfo* pSrc, AMVPInfo* pDst)
{
	pDst->iN = pSrc->iN;
	for (Int i = 0; i < pSrc->iN; i++)
	{
		pDst->m_acMvCand[i] = pSrc->m_acMvCand[i];
	}
}
Void TEncCu::xCopyYuv2Pic(TComPic* rpcPic, UInt uiCUAddr, UInt uiAbsPartIdx, UInt uiDepth, UInt uiSrcDepth)
{
	UInt uiAbsPartIdxInRaster = g_auiZscanToRaster[uiAbsPartIdx];
	UInt uiSrcBlkWidth = rpcPic->getNumPartInCtuWidth() >> (uiSrcDepth);
	UInt uiBlkWidth = rpcPic->getNumPartInCtuWidth() >> (uiDepth);
	UInt uiPartIdxX = ((uiAbsPartIdxInRaster % rpcPic->getNumPartInCtuWidth()) % uiSrcBlkWidth) / uiBlkWidth;
	UInt uiPartIdxY = ((uiAbsPartIdxInRaster / rpcPic->getNumPartInCtuWidth()) % uiSrcBlkWidth) / uiBlkWidth;
	UInt uiPartIdx = uiPartIdxY * (uiSrcBlkWidth / uiBlkWidth) + uiPartIdxX;
	m_ppcRecoYuvBest[uiSrcDepth]->copyToPicYuv(rpcPic->getPicYuvRec(), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);

	m_ppcPredYuvBest[uiSrcDepth]->copyToPicYuv(rpcPic->getPicYuvPred(), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);
}

Void TEncCu::xCopyYuv2Tmp(UInt uiPartUnitIdx, UInt uiNextDepth)
{
	UInt uiCurrDepth = uiNextDepth - 1;
	m_ppcRecoYuvBest[uiNextDepth]->copyToPartYuv(m_ppcRecoYuvTemp[uiCurrDepth], uiPartUnitIdx);
	m_ppcPredYuvBest[uiNextDepth]->copyToPartYuv(m_ppcPredYuvBest[uiCurrDepth], uiPartUnitIdx);
}

/** Function for filling the PCM buffer of a CU using its original sample array
 * \param pCU pointer to current CU
 * \param pOrgYuv pointer to original sample array
 */
Void TEncCu::xFillPCMBuffer(TComDataCU* pCU, TComYuv* pOrgYuv)
{
	const ChromaFormat format = pCU->getPic()->getChromaFormat();
	const UInt numberValidComponents = getNumberValidComponents(format);
	for (UInt componentIndex = 0; componentIndex < numberValidComponents; componentIndex++)
	{
		const ComponentID component = ComponentID(componentIndex);

		const UInt width = pCU->getWidth(0) >> getComponentScaleX(component, format);
		const UInt height = pCU->getHeight(0) >> getComponentScaleY(component, format);

		Pel* source = pOrgYuv->getAddr(component, 0, width);
		Pel* destination = pCU->getPCMSample(component);

		const UInt sourceStride = pOrgYuv->getStride(component);

		for (Int line = 0; line < height; line++)
		{
			for (Int column = 0; column < width; column++)
			{
				destination[column] = source[column];
			}

			source += sourceStride;
			destination += width;
		}
	}
}

#if ADAPTIVE_QP_SELECTION
/** Collect ARL statistics from one block
  */
Int TEncCu::xTuCollectARLStats(TCoeff* rpcCoeff, TCoeff* rpcArlCoeff, Int NumCoeffInCU, Double* cSum, UInt* numSamples)
{
	for (Int n = 0; n < NumCoeffInCU; n++)
	{
		TCoeff u = abs(rpcCoeff[n]);
		TCoeff absc = rpcArlCoeff[n];

		if (u != 0)
		{
			if (u < LEVEL_RANGE)
			{
				cSum[u] += (Double)absc;
				numSamples[u]++;
			}
			else
			{
				cSum[LEVEL_RANGE] += (Double)absc - (Double)(u << ARL_C_PRECISION);
				numSamples[LEVEL_RANGE]++;
			}
		}
	}

	return 0;
}

//! Collect ARL statistics from one CTU
Void TEncCu::xCtuCollectARLStats(TComDataCU* pCtu)
{
	Double cSum[LEVEL_RANGE + 1];     //: the sum of DCT coefficients corresponding to data type and quantization output
	UInt numSamples[LEVEL_RANGE + 1]; //: the number of coefficients corresponding to data type and quantization output

	TCoeff* pCoeffY = pCtu->getCoeff(COMPONENT_Y);
	TCoeff* pArlCoeffY = pCtu->getArlCoeff(COMPONENT_Y);
	const TComSPS& sps = *(pCtu->getSlice()->getSPS());

	const UInt uiMinCUWidth = sps.getMaxCUWidth() >> sps.getMaxTotalCUDepth(); // NOTE: ed - this is not the minimum CU width. It is the square-root of the number of coefficients per part.
	const UInt uiMinNumCoeffInCU = 1 << uiMinCUWidth;                          // NOTE: ed - what is this?

	memset(cSum, 0, sizeof(Double) * (LEVEL_RANGE + 1));
	memset(numSamples, 0, sizeof(UInt) * (LEVEL_RANGE + 1));

	// Collect stats to cSum[][] and numSamples[][]
	for (Int i = 0; i < pCtu->getTotalNumPart(); i++)
	{
		UInt uiTrIdx = pCtu->getTransformIdx(i);

		if (pCtu->isInter(i) && pCtu->getCbf(i, COMPONENT_Y, uiTrIdx))
		{
			xTuCollectARLStats(pCoeffY, pArlCoeffY, uiMinNumCoeffInCU, cSum, numSamples);
		} //Note that only InterY is processed. QP rounding is based on InterY data only.

		pCoeffY += uiMinNumCoeffInCU;
		pArlCoeffY += uiMinNumCoeffInCU;
	}

	for (Int u = 1; u < LEVEL_RANGE; u++)
	{
		m_pcTrQuant->getSliceSumC()[u] += cSum[u];
		m_pcTrQuant->getSliceNSamples()[u] += numSamples[u];
	}
	m_pcTrQuant->getSliceSumC()[LEVEL_RANGE] += cSum[LEVEL_RANGE];
	m_pcTrQuant->getSliceNSamples()[LEVEL_RANGE] += numSamples[LEVEL_RANGE];
}
#endif

#if SVC_EXTENSION
#if N0383_IL_CONSTRAINED_TILE_SETS_SEI
Bool TEncCu::xCheckTileSetConstraint(TComDataCU*& rpcCU)
{
	Bool disableILP = false;

	if (rpcCU->getPic()->getLayerId() == (m_pcEncCfg->getNumLayer() - 1) && m_pcEncCfg->getInterLayerConstrainedTileSetsSEIEnabled() && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(rpcCU->getCtuRsAddr()) >= 0)
	{
		if (rpcCU->getPic()->getPicSym()->getTileSetType(rpcCU->getCtuRsAddr()) == 2)
		{
			disableILP = true;
		}
		if (rpcCU->getPic()->getPicSym()->getTileSetType(rpcCU->getCtuRsAddr()) == 1)
		{
			Int currCUaddr = rpcCU->getCtuRsAddr();
			Int frameWitdhInCU = rpcCU->getPic()->getPicSym()->getFrameWidthInCtus();
			Int frameHeightInCU = rpcCU->getPic()->getPicSym()->getFrameHeightInCtus();
			Bool leftCUExists = (currCUaddr % frameWitdhInCU) > 0;
			Bool aboveCUExists = (currCUaddr / frameWitdhInCU) > 0;
			Bool rightCUExists = (currCUaddr % frameWitdhInCU) < (frameWitdhInCU - 1);
			Bool belowCUExists = (currCUaddr / frameWitdhInCU) < (frameHeightInCU - 1);
			Int currTileSetIdx = rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr);
			// Check if CU is at tile set boundary
			if ((leftCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr - 1) != currTileSetIdx) ||
				(leftCUExists && aboveCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr - frameWitdhInCU - 1) != currTileSetIdx) ||
				(aboveCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr - frameWitdhInCU) != currTileSetIdx) ||
				(aboveCUExists && rightCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr - frameWitdhInCU + 1) != currTileSetIdx) ||
				(rightCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr + 1) != currTileSetIdx) ||
				(rightCUExists && belowCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr + frameWitdhInCU + 1) != currTileSetIdx) ||
				(belowCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr + frameWitdhInCU) != currTileSetIdx) ||
				(belowCUExists && leftCUExists && rpcCU->getPic()->getPicSym()->getTileSetIdxMap(currCUaddr + frameWitdhInCU - 1) != currTileSetIdx))
			{
				disableILP = true; // Disable ILP in tile set boundary CU
			}
		}
	}

	return disableILP;
}

Void TEncCu::xVerifyTileSetConstraint(TComDataCU*& rpcCU)
{
	if (rpcCU->getPic()->getLayerId() == (m_pcEncCfg->getNumLayer() - 1) && m_pcEncCfg->getInterLayerConstrainedTileSetsSEIEnabled() &&
		rpcCU->getPic()->getPicSym()->getTileSetIdxMap(rpcCU->getCtuRsAddr()) >= 0 && m_disableILP)
	{
		UInt numPartitions = rpcCU->getPic()->getNumPartitionsInCtu();
		for (UInt i = 0; i < numPartitions; i++)
		{
			if (!rpcCU->isIntra(i))
			{
				for (UInt refList = 0; refList < 2; refList++)
				{
					if (rpcCU->getInterDir(i) & (1 << refList))
					{
						TComCUMvField* mvField = rpcCU->getCUMvField(RefPicList(refList));
						if (mvField->getRefIdx(i) >= 0)
						{
							assert(!(rpcCU->getSlice()->getRefPic(RefPicList(refList), mvField->getRefIdx(i))->isILR(rpcCU->getPic()->getLayerId())));
						}
					}
				}
			}
		}
	}
}
#endif

#if ENCODER_FAST_MODE
Void TEncCu::xCheckRDCostILRUni(TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt refLayerId)
{
	UChar uhDepth = rpcTempCU->getDepth(0);
	rpcTempCU->setDepthSubParts(uhDepth, 0);
#if SKIP_FLAG
	rpcTempCU->setSkipFlagSubParts(false, 0, uhDepth);
#endif
	rpcTempCU->setPartSizeSubParts(SIZE_2Nx2N, 0, uhDepth); //2Nx2N
	rpcTempCU->setPredModeSubParts(MODE_INTER, 0, uhDepth);
	rpcTempCU->setCUTransquantBypassSubParts(m_pcEncCfg->getCUTransquantBypassFlagForceValue(), 0, uhDepth);
	Bool exitILR = m_pcPredSearch->predInterSearchILRUni(rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth], refLayerId);
	if (!exitILR)
	{
		return;
	}
	m_pcPredSearch->encodeResAndCalcRdInterCU(rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false);
	rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost(rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion());
	xCheckDQP(rpcTempCU);
	xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
	return;
}
#endif

bool TEncCu::isEdge(int depth, int index)
{
	if (depth == 1)
		return index != 0;
	else if (depth == 2)
		return (index % 4 == 3) || (index >= 12);
	return false;
}

double TEncCu::hytest(TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPreYuv, int uiDepth)
{
	if (uiDepth < 1 || uiDepth > 2)
		return -1;
	UInt uiPartSize = pcCU->getWidth(0);
	UInt uiTrUnitIdx = 0, i = 0, j = 0;
	const ComponentID compID = ComponentID(0);
	const Int uiPartWidth = uiPartSize >> pcOrgYuv->getComponentScaleX(compID);
	const Int uiPartHeight = uiPartSize >> pcOrgYuv->getComponentScaleY(compID);
	//int **p=new int*[uiPartHeight];
	//for(int i=0;i<uiPartHeight;i++)             
	//	p[i]=new int[uiPartWidth];
	vector<vector<int>> p(uiPartHeight, vector<int>(uiPartWidth));
	const Pel* pSrc0 = pcOrgYuv->getAddr(compID, uiTrUnitIdx, uiPartWidth);
	const Pel* pSrc1 = pcPreYuv->getAddr(compID, uiTrUnitIdx, uiPartWidth);

	const Int  iSrc0Stride = pcOrgYuv->getStride(compID);
	const Int  iSrc1Stride = pcPreYuv->getStride(compID);
	long n = uiPartWidth * uiPartHeight;

	vector<double> ave_hor(2, 0.0);
	vector<double> ave_ver(2, 0.0);
	for (i = 0; i < uiPartHeight; i++)
	{
		for (j = 0; j < uiPartWidth; j++)
		{
			p[i][j] = abs(pSrc1[j] - pSrc0[j]);
			if (uiDepth < 3)
			{
				if (i < uiPartHeight / 2)
					ave_hor[0] += p[i][j];
				else
					ave_hor[1] += p[i][j];
				if (j < uiPartWidth / 2)
					ave_ver[0] += p[i][j];
				else
					ave_ver[1] += p[i][j];
			}
		}
		pSrc0 += iSrc0Stride;
		pSrc1 += iSrc1Stride;
	}
	ave_ver[0] /= n / 2.0;
	ave_ver[1] /= n / 2.0;
	ave_hor[0] /= n / 2.0;
	ave_hor[1] /= n / 2.0;
	const double value = max(fabs(ave_ver[0] - ave_ver[1]), fabs(ave_hor[0] - ave_hor[1]));
	//vector<double> paras1({ 0.9976, -0.0825, -0.06541, 0.016, -0.001347, 3.888e-05 });
	//vector<double> paras2({ 0.9911, -0.002199, -0.01533, 0.001641, -5.38e-05 });
	double ans = 0;
	if (uiDepth == 1)
	{
		for (int i = 0; i < paras1.size(); i++)
			ans += pow(value, i) * paras1[i];
	}
	else if (uiDepth == 2)
	{
		for (int i = 0; i < paras2.size(); i++)
			ans += pow(value, i) * paras2[i];
	}
	return ans;
}

double TEncCu::getCorrelation(dg& x, dg& y)
{
	double a = 0, b = 0, c = 0;
	double avgx = x.first;
	double avgy = y.first;
	vector<vector<int>>& gx = x.second;
	vector<vector<int>>& gy = y.second;
	int len = x.second.size();
	for (int i = 0; i < len; i++)
	{
		for (int j = 0; j < len; j++)
		{
			a += (gx[i][j] - avgx) * (gy[i][j] - avgy);
			b += (gx[i][j] - avgx) * (gx[i][j] - avgx);
			c += (gy[i][j] - avgy) * (gy[i][j] - avgy);
		}
	}
	return a / (sqrt(b) * sqrt(c));
}

double TEncCu::calD(vector<double>& w, vector<double>& pf, int idx)
{
	double a = 0, b = 0;
	for (int i = 0; i < 3; i++)
	{
		a += w[i] * pf[i];
		b += w[i] * (1 - pf[i]);
	}
	double x = (2 * pf[idx] - 1) * exp(b);
	double y = exp(a) + exp(b);
	return x / y;
}
#endif //SVC_EXTENSION
//! \}
