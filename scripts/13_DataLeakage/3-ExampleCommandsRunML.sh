taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDeg.txt_Balanced_500.txt 1 0 4E_Deg1 &> 4E_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDegD2.txt_Balanced_500.txt 1 0 4E_Deg2 &> 4E_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDeg.txt_Balanced_500.txt 1 0 4F_Deg1 &> 4F_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDegD2.txt_Balanced_500.txt 1 0 4F_Deg2 &> 4F_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDeg.txt_Balanced_500.txt 1 0 4G_Deg1 &> 4G_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDegD2.txt_Balanced_500.txt 1 0 4G_Deg2 &> 4G_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDeg.txt_Balanced_500.txt 1 0 4H_Deg1 &> 4H_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDegD2.txt_Balanced_500.txt 1 0 4H_Deg2 &> 4H_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeq.txt_Balanced_500.txt 1 0 4A_Deg1_W2 &> 4A_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeqD2.txt_Balanced_500.txt 1 0 4A_Deg2_W2 &> 4A_Deg2_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotif.txt_Balanced_500.txt 1 0 4B_Deg1_W2 &> 4B_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotifD2.txt_Balanced_500.txt 1 0 4B_Deg2_W2 &> 4B_Deg2_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgr.txt_Balanced_500.txt 1 0 4C_Deg1_W2 &> 4C_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgrD2.txt_Balanced_500.txt 1 0 4C_Deg2_W2 &> 4C_Deg2_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChem.txt_Balanced_500.txt 1 0 4D_Deg1_W2 &> 4D_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChemD2.txt_Balanced_500.txt 1 0 4D_Deg2_W2 &> 4D_Deg2_W2.txt &




(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDeg.txt_Balanced_500.txt 1 0 4E_Deg1 &> 4E_Deg1.txt &
[1] 3381853
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDegD2.txt_Balanced_500.txt 1 0 4E_Deg2 &> 4E_Deg2.txt &
[2] 3381854
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDeg.txt_Balanced_500.txt 1 0 4F_Deg1 &> 4F_Deg1.txt &
[3] 3381855
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDegD2.txt_Balanced_500.txt 1 0 4F_Deg2 &> 4F_Deg2.txt &
[4] 3381856
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDeg.txt_Balanced_500.txt 1 0 4G_Deg1 &> 4G_Deg1.txt &
[5] 3381857
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDegD2.txt_Balanced_500.txt 1 0 4G_Deg2 &> 4G_Deg2.txt &
[6] 3381858
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDeg.txt_Balanced_500.txt 1 0 4H_Deg1 &> 4H_Deg1.txt &
[7] 3381859
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDegD2.txt_Balanced_500.txt 1 0 4H_Deg2 &> 4H_Deg2.txt &
[8] 3381860
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeq.txt_Balanced_500.txt 1 0
4A_Deg1_W2 &> 4A_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeqD2.txt_Balanced_500.txt 1 0 4A_Deg2_W2 &> 4A_Deg2_W2.txt &
[9] 3381861
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeqD2.txt_Balanced_500.txt 1
0 4A_Deg2_W2 &> 4A_Deg2_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotif.txt_Balanced_500.txt 1 0 4B_Deg1_W2 &> 4B_Deg1_W2.txt &
[10] 3381862
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotif.txt_Balanced_500.txt
1 0 4B_Deg1_W2 &> 4B_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotifD2.txt_Balanced_500.txt 1 0 4B_Deg2_W2 &> 4B_Deg2_W2.txt &
[11] 3381863
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotifD2.txt_Balanced_500.txt 1 0 4B_Deg2_W2 &> 4B_Deg2_W2.txt &
[12] 3381864
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgr.txt_Balanced_500.txt 1 0
4C_Deg1_W2 &> 4C_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgrD2.txt_Balanced_500.txt 1 0 4C_Deg2_W2 &> 4C_Deg2_W2.txt &
[13] 3381866
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgrD2.txt_Balanced_500.txt 1
0 4C_Deg2_W2 &> 4C_Deg2_W2.txt &
[14] 3381867
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChem.txt_Balanced_500.txt 1
0 4D_Deg1_W2 &> 4D_Deg1_W2.txt &
[15] 3381868
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4Mem.py 6 500 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChemD2.txt_Balanced_500.txt
1 0 4D_Deg2_W2 &> 4D_Deg2_W2.txt &
[16] 3381869









taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDeg.txt_Balanced_500.txt 1 0 4TRE_Deg1 &> Transformer4TRE_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDegD2.txt_Balanced_500.txt 1 0 4TRE_Deg2 &> Transformer4TRE_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDeg.txt_Balanced_500.txt 1 0 4TRF_Deg1 &> Transformer4TRF_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDegD2.txt_Balanced_500.txt 1 0 4TRF_Deg2 &> Transformer4TRF_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDeg.txt_Balanced_500.txt 1 0 4TRG_Deg1 &> Transformer4TRG_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDegD2.txt_Balanced_500.txt 1 0 4TRG_Deg2 &> Transformer4TRG_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDeg.txt_Balanced_500.txt 1 0 4TRH_Deg1 &> Transformer4TRH_Deg1.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDegD2.txt_Balanced_500.txt 1 0 4TRH_Deg2 &> Transformer4TRH_Deg2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeq.txt_Balanced_500.txt 1 0 4TRA_Deg1_W2 &> Transformer4TRA_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeqD2.txt_Balanced_500.txt 1 0 4TRA_Deg2_W2 &> Transformer4TRA_Deg2_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotif.txt_Balanced_500.txt 1 0 4TRB_Deg1_W2 &> Transformer4TRB_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotifD2.txt_Balanced_500.txt 1 0 4TRB_Deg2_W2 &> Transformer4TRB_Deg2_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgr.txt_Balanced_500.txt 1 0 4TRC_Deg1_W2 &> Transformer4TRC_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgrD2.txt_Balanced_500.txt 1 0 4TRC_Deg2_W2 &> Transformer4TRC_Deg2_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChem.txt_Balanced_500.txt 1 0 4TRD_Deg1_W2 &> Transformer4TRD_Deg1_W2.txt &
taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChemD2.txt_Balanced_500.txt 1 0 4TRD_Deg2_W2 &> Transformer4TRD_Deg2_W2.txt &


(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDeg.txt_Balanced_500.txt 1 0 4TRE_Deg1 &> Transformer4TRE_Deg1.txt &
[1] 3112088
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_E_EpiSeq_ParaSeq_NoDegD2.txt_Balanced_500.txt 1 0 4TRE_Deg2 &> Transformer4TRE_Deg2.txt &
[2] 3112089
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDeg.txt_Balanced_500.txt 1 0 4TRF_Deg1 &> Transformer4TRF_Deg1.txt &
[3] 3112090
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_F_EpiMotif_ParaMotif_NoDegD2.txt_Balanced_500.txt 1 0 4TRF_Deg2 &> Transformer4TRF_Deg2.txt &
[4] 3112091
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDeg.txt_Balanced_500.txt 1 0 4TRG_Deg1 &> Transformer4TRG_Deg1.txt &
[5] 3112092
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_G_EpiAgr_ParaAgr_NoDegD2.txt_Balanced_500.txt 1 0 4TRG_Deg2 &> Transformer4TRG_Deg2.txt &
[6] 3112093
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDeg.txt_Balanced_500.txt 1 0 4TRH_Deg1 &> Transformer4TRH_Deg1.txt &
[7] 3112094
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 1 1 ../Task4/Task4_H_EpiChem_ParaChem_NoDegD2.txt_Balanced_500.txt 1 0 4TRH_Deg2 &> Transformer4TRH_Deg2.txt &
[8] 3112095
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeq.txt_Balanced_500.txt 1 0 4TRA_Deg1_W2 &> Transformer4TRA_Deg1_W2.txt &
[9] 3112096
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_A_EpiSeq_ParaSeqD2.txt_Balanced_500.txt 1 0 4TRA_Deg2_W2 &> Transformer4TRA_Deg2_W2.txt &
[10] 3112097
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotif.txt_Balanced_500.txt 1 0 4TRB_Deg1_W2 &> Transformer4TRB_Deg1_W2.txt &
[11] 3112098
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_B_EpiMotif_ParaMotifD2.txt_Balanced_500.txt 1 0 4TRB_Deg2_W2 &> Transformer4TRB_Deg2_W2.txt &
[12] 3112099
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgr.txt_Balanced_500.txt 1 0 4TRC_Deg1_W2 &> Transformer4TRC_Deg1_W2.txt &
[13] 3112100
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_C_EpiAgr_ParaAgrD2.txt_Balanced_500.txt 1 0 4TRC_Deg2_W2 &> Transformer4TRC_Deg2_W2.txt &
[14] 3112101
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChem.txt_Balanced_500.txt 1 0 4TRD_Deg1_W2 &> Transformer4TRD_Deg1_W2.txt &
[15] 3112102
(tensorflow2GPU) [pprobert@immunohub01 PyMem]$ taskset --cpu-list 150-220 nice -n 19 ./Task4-TransformerMem.py 16 512 10000 10 1 2 2 ../Task4/Task4_D_EpiChem_ParaChemD2.txt_Balanced_500.txt 1 0 4TRD_Deg2_W2 &> Transformer4TRD_Deg2_W2.txt &
[16] 3112103
