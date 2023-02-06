## scanpy

[https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.obs_df.html#scanpy.get.obs_df](https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.obs_df.html#scanpy.get.obs_df)

- Load a scanpy object:

```
adata = sc.read_h5ad('anndata.h5ad')
```

```
col_ids = ["AAACCCACACTCCTGT_TSP6_Liver_NA_10X_1_1", "AAACGAAGTACCAGAG_TSP6_Liver_NA_10X_1_1"] test_subset = sc.get.var_df(adata, keys=col_ids)
```

- Get gene names:

```
adata.var_names
```

- Get cell names:

```
print(set(adata.obs_names))
```

- Number of cells:

```
adata.shape[0]
```

- Get metadata:

```
adata_original.obs
```

- Get raw data:
```
adata2 = adata.raw.to_adata()
```

- Loading large dataset:
  + large datasets in h5ad format can be more easily handled if we read them using and additional ” backed= ‘r’ ” to the anndata.read function. 
  + We then can make subsets of cells (genes) with only the .var and .obs in memory, and keep the data itself stored on disk
  + reference: https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_h5ad.html#scanpy-read-h5ad
