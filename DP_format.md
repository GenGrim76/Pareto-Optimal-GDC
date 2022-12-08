## DP output files Documentation

**PLEASE NOTE:** Assume `D` the genomic input dataset. We assume `DSK_Dict` the Dictionary obtained from DSK in the (`di` `fi`) format, where `di` is a canonical k-mer present in D, and `fi` is the relative frequency. DSK We assume `Dk` to be the dictionary of k-mers obtained from DSK, and `Fk` the dictionary of relative frequencies. All data related to compression are expressed in [bytes], while those related to post-processing time are expressed in [ms].

***Base Case Scenario***

```
DP_Base_size.txt
uncompressed_D, bzip2(D), lz4(D), zstd(D), MFC(D), SPRING(D)
uncompressed_D, FM-index(D)

DP_Base_time.txt
t_compr_bzip2(D), t_compr_lz4(D), t_compr_zstd(D), t_compr_MFC(D), t_compr_SPRING(D)
t_decompr_bzip2(D), t_decompr_lz4(D), t_decompr_zstd(D), t_decompr_MFC(D), t_decompr_SPRING(D)
t_compr_FM-index(D)
t_decompr_FM-index(D)
```


***SD-RAM Scenario***




***CD-NRAM Scenario***

```
DP[0/1/2]_k[x]-size.txt
uncompressed_Dk, bzip2(Dk), lz4(Dk), zstd(Dk), MFC(Dk), SPRING(Dk) - uncompressed_Dk, FM-index(Dk)
uncompressed_Fk, bzip2(Fk), lz4(Fk), zstd(Fk), BIC(Fk), Opt-PFOR(Fk) - uncompressed(DSK_Dict), BCSF(DSK_Dict), bzip2(BCSF(DSK_Dict)), lz4(BCSF(DSK_Dict)), zstd(BCSF(DSK_Dict))
uncompressed_Gk, bzip2(Gk), lz4(Gk), zstd(Gk), BIC(Gk), Opt-PFOR(Gk) 

DP[0/1/2]_k[x]-time.txt
DSK
DSK to {Dk, Fk}
t_compr_bzip2(Dk), t_compr_lz4(Dk), t_compr_zstd(Dk), t_compr_MFC(Dk), t_compr_SPRING(Dk) - t_compr_FM-index(Dk)
t_decompr_bzip2(Dk), t_decompr_lz4(Dk), t_decompr_zstd(Dk), t_decompr_MFC(Dk), t_decompr_SPRING(Dk) - t_decompr_FM-index(Dk)
t_compr_bzip2(Fk), t_compr_lz4(Fk), t_compr_zstd(Fk), t_compr_BIC(Fk), t_compr_Opt-PFOR(Fk) - t_build_BCSF(DSK_Dict), bzip2(BCSF(DSK_Dict)), lz4(BCSF(DSK_Dict)), zstd(BCSF(DSK_Dict))
t_decompr_bzip2(Fk), t_decompr_lz4(Fk), t_decompr_zstd(Fk), t_decompr_BIC(Fk), t_decompr_Opt-PFOR(Fk) - t_decompr_bzip2(BCSF(DSK_Dict)), t_decompr_lz4(BCSF(DSK_Dict)), t_decompr_zstd(BCSF(DSK_Dict))
Fk to Gk
t_compr_bzip2(Gk), t_compr_lz4(Gk), t_compr_zstd(Gk), t_compr_BIC(Gk), t_compr_Opt-PFOR(Gk)
t_decompr_bzip2(Gk), t_decompr_lz4(Gk), t_decompr_zstd(Gk), t_decompr_BIC(Gk), t_decompr_Opt-PFOR(Gk)
Gk to Fk
```
