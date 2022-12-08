## DP output files Documentation

**PLEASE NOTE:** 

Assume the following: 
- `D` is the genomic input dataset 
- `DSK_Dict` is the dictionary obtained from DSK in the (`di` `fi`) format, where `di` is a canonical k-mer present in D, and `fi` is the relative frequency
- `Dk` is the dictionary of k-mers obtained from DSK
- `Fk` is the dictionary of relative frequencies of 'di' in 'Dk'.
- `S` is the Spectrum Preserving String Set obtained from ESScompress. 
 
All data related to compression are expressed in [bytes], while those related to post-processing time are expressed in [ms]. For the convenience of the reader, we indicate the progressive line number.


***Base Case Scenario***

```
DP_Base_size.txt
1: uncompressed_D, bzip2(D), lz4(D), zstd(D), MFC(D), SPRING(D)
2: uncompressed_D, FM-index(D)

DP_Base_time.txt
1: t_compr_bzip2(D), t_compr_lz4(D), t_compr_zstd(D), t_compr_MFC(D), t_compr_SPRING(D)
2: t_decompr_bzip2(D), t_decompr_lz4(D), t_decompr_zstd(D), t_decompr_MFC(D), t_decompr_SPRING(D)
3: t_compr_FM-index(D)
4: t_decompr_FM-index(D)
```

When we are in the exact solution case, we consider the following:

***SD-RAM Scenario***

```
DP3_k[x]-size-SD_RAM.txt
1: uncompressed_S, FM-index(S)


DP3_k[x]-time-SD_RAM.txt
1: t_compr_FM-index(S)
2: t_decompr_FM-index(S)
3: t_recovery(Dk) from S
4: t_prepare(Dk) for recovery Fk
5: t_recovery(Fk) from Dk and BCSF
```



***CD-NRAM Scenario***

```
DP0_k[x]-size.txt
1: uncompressed_Dk, bzip2(Dk), lz4(Dk), zstd(Dk), MFC(Dk), SPRING(Dk) - uncompressed_Dk, FM-index(Dk)
2: uncompressed_Fk, bzip2(Fk), lz4(Fk), zstd(Fk), BIC(Fk), Opt-PFOR(Fk) - uncompressed(DSK_Dict), BCSF(DSK_Dict), bzip2(BCSF(DSK_Dict)), lz4(BCSF(DSK_Dict)), zstd(BCSF(DSK_Dict))
3: uncompressed_Gk, bzip2(Gk), lz4(Gk), zstd(Gk), BIC(Gk), Opt-PFOR(Gk) 

DP0_k[x]-time.txt
1: DSK
2: DSK to {Dk, Fk}
3: t_compr_bzip2(Dk), t_compr_lz4(Dk), t_compr_zstd(Dk), t_compr_MFC(Dk), t_compr_SPRING(Dk) - t_compr_FM-index(Dk)
4: t_decompr_bzip2(Dk), t_decompr_lz4(Dk), t_decompr_zstd(Dk), t_decompr_MFC(Dk), t_decompr_SPRING(Dk) - t_decompr_FM-index(Dk)
5: t_compr_bzip2(Fk), t_compr_lz4(Fk), t_compr_zstd(Fk), t_compr_BIC(Fk), t_compr_Opt-PFOR(Fk) - t_build_BCSF(DSK_Dict), bzip2(BCSF(DSK_Dict)), lz4(BCSF(DSK_Dict)), zstd(BCSF(DSK_Dict))
6: t_decompr_bzip2(Fk), t_decompr_lz4(Fk), t_decompr_zstd(Fk), t_decompr_BIC(Fk), t_decompr_Opt-PFOR(Fk) - t_decompr_bzip2(BCSF(DSK_Dict)), t_decompr_lz4(BCSF(DSK_Dict)), t_decompr_zstd(BCSF(DSK_Dict))
7: Fk to Gk
8: t_compr_bzip2(Gk), t_compr_lz4(Gk), t_compr_zstd(Gk), t_compr_BIC(Gk), t_compr_Opt-PFOR(Gk)
9: t_decompr_bzip2(Gk), t_decompr_lz4(Gk), t_decompr_zstd(Gk), t_decompr_BIC(Gk), t_decompr_Opt-PFOR(Gk)
10: Gk to Fk

DP[1/2]_k[x]-size.txt
1: uncompressed_Dk, bzip2(Dk), lz4(Dk), zstd(Dk), MFC(Dk), SPRING(Dk) - uncompressed_Dk, FM-index(Dk)
2: uncompressed_Fk, bzip2(Fk), lz4(Fk), zstd(Fk), BIC(Fk), Opt-PFOR(Fk)
3: uncompressed_Gk, bzip2(Gk), lz4(Gk), zstd(Gk), BIC(Gk), Opt-PFOR(Gk) 

DP[1/2]_k[x]-time.txt
1: DSK
2: DSK to {Dk, Fk}
3: t_compr_bzip2(Dk), t_compr_lz4(Dk), t_compr_zstd(Dk), t_compr_MFC(Dk), t_compr_SPRING(Dk) - t_compr_FM-index(Dk)
4: t_decompr_bzip2(Dk), t_decompr_lz4(Dk), t_decompr_zstd(Dk), t_decompr_MFC(Dk), t_decompr_SPRING(Dk) - t_decompr_FM-index(Dk)
5: t_compr_bzip2(Fk), t_compr_lz4(Fk), t_compr_zstd(Fk), t_compr_BIC(Fk), t_compr_Opt-PFOR(Fk)
6: t_decompr_bzip2(Fk), t_decompr_lz4(Fk), t_decompr_zstd(Fk), t_decompr_BIC(Fk), t_decompr_Opt-PFOR(Fk)
7: Fk to Gk
8: t_compr_bzip2(Gk), t_compr_lz4(Gk), t_compr_zstd(Gk), t_compr_BIC(Gk), t_compr_Opt-PFOR(Gk)
9: t_decompr_bzip2(Gk), t_decompr_lz4(Gk), t_decompr_zstd(Gk), t_decompr_BIC(Gk), t_decompr_Opt-PFOR(Gk)
10: Gk to Fk


DP3_k[x]-size-CD_NRAM.txt
1: uncompressed_Fk, bzip2(Fk), lz4(Fk), zstd(Fk), BIC(Fk), Opt-PFOR(Fk)
2: uncompressed_Gk, bzip2(Gk), lz4(Gk), zstd(Gk), BIC(Gk), Opt-PFOR(Gk)
3: uncompressed_S, bzip2(S), lz4(S), zstd(S), MFC(S), SPRING(S)


DP3_k[x]-time-CD_NRAM.txt
1: DSK
2: ESSCompress
3: F'k Bijection
4: t_compr_bzip2(Fk), t_compr_lz4(Fk), t_compr_zstd(Fk), t_compr_BIC(Fk), t_compr_Opt-PFOR(Fk)
5: t_decompr_bzip2(Fk), t_decompr_lz4(Fk), t_decompr_zstd(Fk), t_decompr_BIC(Fk), t_decompr_Opt-PFOR(Fk) 
6: Fk to Gk
7: t_compr_bzip2(Gk), t_compr_lz4(Gk), t_compr_zstd(Gk), t_compr_BIC(Gk), t_compr_Opt-PFOR(Gk)
8: t_decompr_bzip2(Gk), t_decompr_lz4(Gk), t_decompr_zstd(Gk), t_decompr_BIC(Gk), t_decompr_Opt-PFOR(Gk) 
9: Gk to Fk
10: D'k Bijection  
11: t_compr_bzip2(S), t_compr_lz4(S), t_compr_zstd(S), t_compr_MFC(S), t_compr_SPRING(S)
12: t_decompr_bzip2(S), t_decompr_lz4(S), t_decompr_zstd(S), t_decompr_MFC(S), t_decompr_SPRING(S)

```

When we are in the approximate solution case, we consider, together with the above, the following DP3 subcases:

***CD-NRAM Scenario***
```
DP3_k[x]-time-CD_NRAM.txt
1: Approximate Solution
2: t_compr_bzip2(Fk), t_compr_lz4(Fk), t_compr_zstd(Fk), t_compr_BIC(Fk), t_compr_Opt-PFOR(Fk)
3: t_decompr_bzip2(Fk), t_decompr_lz4(Fk), t_decompr_zstd(Fk), t_decompr_BIC(Fk), t_decompr_Opt-PFOR(Fk) 
4: Fk to Gk
5: t_compr_bzip2(Gk), t_compr_lz4(Gk), t_compr_zstd(Gk), t_compr_BIC(Gk), t_compr_Opt-PFOR(Gk)
6: t_decompr_bzip2(Gk), t_decompr_lz4(Gk), t_decompr_zstd(Gk), t_decompr_BIC(Gk), t_decompr_Opt-PFOR(Gk) 
7: Gk to Fk
8: t_compr_bzip2(S), t_compr_lz4(S), tcompr_zstd(S), t_compr_MFC(S), t_compr_SPRING(S)
9: t_decompr_bzip2(S), t_decompr_lz4(S), t_decompr_zstd(S), t_decompr_MFC(S), t_decompr_SPRING(S)
10: DSK
11: ESSCompress
12: F'k Bijection
13: D'k Bijection
14: Bijection Create Approximate Dictionary (Dk, Fk)
```

***SD-RAM Scenario***
```
DP3_k[x]-time-SD_RAM.txt
1: Approximate Solution
2: t_compr_FM-index(S)
3: t_decompr_FM-index(S)
3: t_recovery(Dk) from S
4: t_prepare(Dk) for recovery Fk
5: t_recovery(Fk) from Dk and BCSF
```
