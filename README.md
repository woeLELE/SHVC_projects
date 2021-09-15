# SHVC_projects
该项目是在导师的指导下完成的，用来实现导师提出的视频编码算法，具体的C++代码是我独立完成的。简单的说，就是我帮导师做项目，导师自己写论文，论文题目是《A New Fast Intra Prediction Algorithm for Spatial SHVC》，已送审IEEE。该项目的工作可以分为三个部分：基于残差图像的两部分的显著性差异判断是否将ILR模式作为最优模式，以跳过进行帧内预测的步骤；对ILR模式残生的RDCost进行GMM变换以判断是否将ILR模式作为最优模式，同样为了跳过帧内预测的步骤；如果前两步都不能确定将ILR作为最优模式，那么进行帧内预测，并对帧内预测过程进行优化以降低运算复杂度。

因为编码器本身的工程文件非常大，所以这里我还是说一下我的代码的实现位置

GMM类是我自己实现的，用来根据当前CU的ILR模式的RDCost以及相邻CU的最优RDCost进行GMM变换判断当前CU更可能属于ILR模式还是帧内模式

TEncCu类是编码器原本就存在的类，我做的工作主要从第873开始，并且在xCheckRDCostILRUni方法中做了修改，调用hytest方法以得到split的结果，hytest方法是我在该类中添加的方法，用来计算残差图像的显著性差异。

紧接着后面就是GMM变换的工作，根据ILRCost和相邻CU的RDCost判断当前CU的预测模式。

帧内预测方法是xCheckRDCostIntra，在1143行，在该方法中在1979行调用了estIntraPredLumaQT方法，用来进行亮度帧内预测，该方法的定义在TEncSearch类文件中，而我对亮度的帧内预测过程的优化集中于TEncSearch类文件的2314行后面
