# SHVC_projects
该项目是在导师的指导下完成的，用来实现导师和本人提出的视频编码算法，具体的C++代码是我独立完成的。根据此项目，完成论文《Efficient Hybrid Strategies for Intra Prediction of Spatial SHVC》，已送审IEEE。
该项目的工作可以分为三个部分：对ILR模式的RDCost进行GMM变换以判断是否将ILR模式作为最优模式，如果是就跳过帧内预测；如果不能确定将ILR作为最优模式，那么进行帧内预测，并对帧内预测过程进行优化以降低运算复杂度。随后，如果当前编码单元（CU）的深度小于3，则结合最优编码模式的率失真代价、残差系数以及当前CU与相邻CU之间的相关性程度，通过条件随机场迭代出当前CU终止划分的概率，从而跳过编码子CU。经过基于标准测试环境的实验，编码速度提升接近 75%。
本人的工作集中在TencCu.cpp中的compressCtu（339行）和xCompressCU（593）方法，TEncSearch.cpp中的estIntraPredLumaQT（2187行）方法，以及GMM类和ITEW类的实现等
