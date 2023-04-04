# cellphonedb

## 01 文件准备

* 刘老师的数据，大群数据是相比亚群数目要多了很多，因为去除了部分双细胞，所以在进行文件准备的时候，需要单独多每一个亚群进行提取后再合并。
* Myeloids

```R
library(Seurat)
library(tidyverse)
library(dplyr)

setwd('/root/wangje/Project/刘老师/CPDB/')
load('../Myeloids/Data/Myeloids_CCA.RData')

# count 信息
counts <- as.matrix(scRNA_seurat@assays$RNA@data)
# write.table(counts, file = 'cpdb_Myeloids_counts.txt', sep ='\t', quote=F)

fwrite(counts, file = '/root/wangje/Project/刘老师/CPDB/cpdb_Myeloids_counts.txt',row.names=F)
# meta
meta <- scRNA_seurat@meta.data %>% select(celltype) %>% rownames_to_column('cell')
write.table(meta, file = '/root/wangje/Project/刘老师/CPDB/cpdb_Myeloids_meta.txt',sep ='\t', row.names=F, quote=F)



```



* NK & T cells

```R
library(Seurat)
library(tidyverse)
library(dplyr)
library(vroom)
library(data.table)
load('/root/wangje/Project/刘老师/NK_T/Data/new_NKT_CCA.RData')
# count
counts <- as.matrix(scRNA_seurat@assays$RNA@data)

# write.table(counts, file = '/root/wangje/Project/刘老师/CPDB/cpdb_NKT_counts.txt', sep ='\t', quote=F)
fwrite(counts, file = '/root/wangje/Project/刘老师/CPDB/cpdb_NKT_counts.txt',row.names=F, sep ='\t', quote=F)

# vroom
vroom::vroom_write(counts, file = '/root/wangje/Project/刘老师/CPDB/cpdb_NKT_counts_2.txt', col_names = T)
# meta
meta <- scRNA_seurat@meta.data %>% select(celltype) %>% rownames_to_column('cell')
write.table(meta, file = '/root/wangje/Project/刘老师/CPDB/cpdb_NKT_meta.txt',sep ='\t', row.names=F, quote=F)

```

* Fib and Endo

```R
library(Seurat)
library(tidyverse)
library(data.table)
# load
load('/root/wangje/Project/刘老师/合并Endothelials和Fibroblasts/CCA/Data/new_scRNA_CCA.RData')
# counts
counts <- as.matrix(scRNA_seurat@assays$RNA@data)
# write
fwrite(counts, file = '/root/wangje/Project/刘老师/CPDB/cpdb_stroma_counts.txt',sep ='\t', row.names=F)

# meta
meta <- scRNA_seurat@meta.data %>% select(celltype) %>% rownames_to_column('cell')
write.table(meta, file = '/root/wangje/Project/刘老师/CPDB/cpdb_stroma_meta.txt',sep ='\t', row.names=F, quote=F)

```

* Epithelials

```R
library(Seurat)
library(tidyverse)
library(data.table)
#load
load("/root/wangje/Project/刘老师/EPi/新结果/添加肝脏/new/new_FastMNN.RData")
# count
counts <- as.matrix(scRNA_seurat@assays$RNA@data)
# write
fwrite(counts, file = '/root/wangje/Project/刘老师/CPDB/cpdb_Epithelials_counts.txt',sep ='\t', row.names=F)

# meta
scRNA_seurat@meta.data[which(scRNA_seurat$RNA_snn_res.0.05 %in% 1), 'celltype'] <- 'Proliferating cells'
scRNA_seurat@meta.data[which(scRNA_seurat$RNA_snn_res.0.05 %in% 0), 'celltype'] <- 'cancer cells'
meta <- scRNA_seurat@meta.data %>% select(celltype) %>% rownames_to_column('cell')
write.table(meta, file = '/root/wangje/Project/刘老师/CPDB/cpdb_Epithelials_meta.txt',sep ='\t', row.names=F, quote=F)

```

* B cells

```r
library(Seurat)
library(tidyverse)
library(dplyr)
library(data.table)

load('/root/wangje/Project/刘老师/Bcells/Data/new_Bcells.RData')
scRNA_seurat <- Bcells
# count
counts <- as.matrix(scRNA_seurat@assays$RNA@data)
# write
fwrite(counts, file = '/root/wangje/Project/刘老师/CPDB/cpdb_Bcells_counts.txt',sep ='\t', row.names=F)
# meta
meta <- scRNA_seurat@meta.data %>% select(celltype) %>% rownames_to_column('cell')
write.table(meta, file = '/root/wangje/Project/刘老师/CPDB/cpdb_Bcells_meta.txt',sep ='\t', row.names=F, quote=F)

# 写出gene的文件
 write.table(rownames(Bcells@assays$RNA@data), file = '/root/wangje/Project/刘老师/CPDB/genes.txt', quote=F, sep ='\t', row.names=F)

```

* 合并文件

#### CellPhoneDB使用4.0版本的时候，运行脚本会报错，显示无法找到数据库，在文件夹中已经下载了数据库。

```shell
paste -d $'\t' genes.txt  cpdb_NKT_counts.txt cpdb_Myeloids_counts.txt cpdb_Bcells_counts.txt cpdb_stroma_counts.txt cpdb_Epithelials_counts.txt > countss.txt
```

```shell
cellphonedb method statistical_analysis meta.txt counts.txt --counts-data=gene_name --threads=40 #这里使用的是cellphonedb3.0的版本
```

* --output-path  Directory where the results will be allocated (the directory must exist) [out]





## 将cellphonedb按照样本进行分析，减少运行时间

### 准备文件

```R
library(data.table)
library(dplyr)
library(tidyverse)
library(future)
library(foreach)
library(doParallel)

counts <- fread('/root/wangje/Project/刘老师/CPDB/counts.txt', header=T,sep='\t', nThread = 90)
info = read.table('/root/wangje/Project/刘老师/CPDB/meta.txt', header=T , sep = '\t')
dup <- info[duplicated(info$cell),1:2]
counts_new <- counts %>% select(which(!colnames(counts) %in% dup$cell))

# 上面的程序写成并行的程序
# cl <- makePSOCKcluster(12)
# registerDoParallel(cl)

# foreach(sp = as.vector(sample$Patient)) %dopar% {
#     tmp <- counts[1:nrow(counts), c(1, grep(paste0('^', sp, '_[A|T|G|C].*'), colnames(counts)))]
#     out_file <- paste0(sp, '_counts.txt')
#     fwrite(tmp, file = out_file, sep = '\t', quote = F)
# }

# stopCluster(cl)
####################################################################################################
###########  筛选每个样本的count信息
sample = read.table('/root/wangje/Project/刘老师/CPDB/sample/sample_ID.txt', header=T, sep = '\t')
for(sp in as.vector(sample$Patient)){
    print(sp)
    tmp = counts %>% select(1,grep(paste0('^', sp, '_[A|T|G|C].*'), colnames(counts)))
    out_file <- paste0(sp, '_counts.txt')
    fwrite(tmp, file = out_file, sep = '\t', quote = F)
}
# 查看列数
ls *.txt | while read id; do awk -F '\t' '{print NF; exit}' $id;done

########### 筛选每个样本的meta信息
meta <- read.table('../meta.txt', header=T,sep= '\t')
meta_new <- meta %>% filter(meta$cell %in% meta$cell[which(!meta$cell %in% dup$cell)])
meta_new$sample <- str_split_fixed(meta_new$cell,'_[A|T|C|G].*', n= 2)[,1]
for(i in unique(meta_new$sample)){
    outfile = paste0(i,'_meta.txt')
    tmp = meta_new %>% filter(sample == i) %>% select(1,2)
    print(head(tmp))
    write.table(tmp, file = outfile, sep ='\t', row.names=F, quote=F)
}
```

### 运行cellphondb

```bash
########### 运行cellphonedb
Pfifo="/tmp/$$.fifo"   # 以PID为名, 防止创建命名管道时与已有文件重名，从而失败
mkfifo $Pfifo          # 创建命名管道
exec 15<>$Pfifo         # 以读写方式打开命名管道, 文件标识符fd为6
                       # fd可取除0, 1, 2,5外0-9中的任意数字
rm -f $Pfifo           # 删除文件, 也可不删除, 不影响后面操作

for i in $(ls -d */)
do
{
    echo ${i%%/*}
    cellphonedb method statistical_analysis ${i %%/*}_meta.txt ${i %%/*}_counts.txt --counts-data=gene_name --threads=40  --output-path ${i} --iterations 10000
}&
done >&15
wait # 等待所有任务结束
exec 15>&- # 关闭管道
```



