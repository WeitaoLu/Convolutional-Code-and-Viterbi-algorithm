# (2,1,2)卷积码实验
# 代码清单
# cc.m 主函数，可直接执行
1. mycc2_1_2(input_data,D),对input_data用卷积码编码,D 为寄存器
2. viterbi_hard(x) 硬判决的viterbi算法
3. viterbi_short_hard(x,L) 深度L的viterbi截短译码算法
4. adam_viterbi(x) 自适应viterbi截短译码算法
# QPSK_Modulation.m 实验一中函数改进版本，包含调制解调和瑞利加性白噪声信道