## DP output files Documentation

**PLEASE NOTE:** Assume D the genomic input dataset.All data related to compression are expressed in [bytes], while those related to post-processing time are expressed in [ms].

***Base Case Scenario***

```
##DP_Base_size.txt##
uncompressed_D, bzip2(D), lz4(D), zstd(D), MFC(D), SPRING(D)
uncompressed_D, FM-index(D)

##DP_Base_time.txt##
t_compr_bzip2(D), t_compr_lz4(D), t_compr_zstd(D), t_compr_MFC(D), t_compr_SPRING(D)
t_decompr_bzip2(D), t_decompr_lz4(D), t_decompr_zstd(D), t_decompr_MFC(D), t_decompr_SPRING(D)
t_compr_FM-index(D)
t_decompr_FM-index(D)
```


***SD-RAM Scenario***




***CD-NRAM Scenario***
