ctMRS论文表格文件说明


一、正文主表

Table 1.csv
名称建议：Study sample information
记录数：421
用途：提供纳入突变分析样本的基础信息，一行代表一个样本。
主要字段：sample_id为样本编号；alias为样本别名；subject_id为受试者编号；biological_sex为生物学性别；phenotype为表型分类。
使用位置：Methods中的研究队列说明，以及正文或补充材料中的样本组成描述。
注意：phenotype包含肺癌、其他癌症、良性疾病和无基线癌症等多个类别，不能直接把421例全部描述为肺癌患者；正式汇总前应按照论文纳入标准重新分类。


Table 2.csv
名称建议：Cell-type-specific ctMRS regulatory module results
记录数：2,932
用途：汇总每个“转录因子–细胞类型”调控模块的CTS、MGRS和MRS/ctMRS分析结果，一行代表一个TF–Celltype组合。
主要字段：TF和Celltype为转录因子及细胞类型；N_targets为靶基因数；CTS_raw/CTS_z为细胞类型特异性原始及标准化分数；N_mutation_linked_targets和MutationLinkedFraction为突变关联靶基因数及比例；MGRS_raw/MGRS_z为突变相关调控评分；Penalty和MRS_raw为惩罚项及最终模块评分；CTS_P、MGRS_P、MRS_P及对应FDR为统计检验结果；Significant、HighConfidence、Candidate和StrongCandidate为筛选层级标记。
使用位置：正文候选调控模块的核心结果表，也可用于生成TF–细胞类型排名图和筛选E2F3等候选模块。
注意：论文中使用Candidate或StrongCandidate时，应明确筛选规则，不能把这些标记自动表述为FDR显著。


Table 3.csv
名称建议：Six-gene LASSO Cox model coefficients
记录数：6
用途：展示六基因预后模型中每个基因的LASSO-Cox回归系数及风险方向。
主要字段：cohort为TCGA-LUAD；endpoint为OS；method为LASSO Cox regression；lambda_rule为lambda.min；rank_by_absolute_coefficient为按系数绝对值排序的名次；gene为模型基因；coefficient为风险评分系数；direction区分Risk和Protective。
使用位置：正文预后模型构建部分，可用于说明CPS1、NPAS2、PTPRH、IGF2BP3和KIAA0319为正向风险基因，COL4A4为保护性基因。


Table 4.csv
名称建议：Prognostic validation of the six-gene signature
记录数：7
用途：汇总六基因风险评分在TCGA-LUAD发现队列和独立GEO验证队列中的总体生存分析结果，一行代表一个队列。
主要字段：dataset和source_type为数据集及发现/验证属性；endpoint为OS；n和events为样本数及死亡事件数；low_risk_n和high_risk_n为两组样本数；HR_high_vs_low、lower95和upper95为HR及95%置信区间；cox_p和logrank_p为两类P值；median_risk_score_cutoff为中位数截点；FDR_logrank为BH校正结果；genes为六基因组成。
使用位置：正文Figure 6对应的预后验证汇总表，以及Results中各独立队列HR、95% CI和P值的依据。
注意：正文报告P值时应明确使用cox_p还是logrank_p，并在全文保持一致。


二、补充表
----------

Supplementary Table 1.csv
名称建议：Gene-level mutation summary
记录数：18,453
用途：将原始逐变异记录按基因汇总，一行代表一个发生突变的基因，替代体积过大的36.8万行逐变异表。
主要字段：gene为基因符号；mutated_sample_count为突变样本数；total_cohort_samples为分母421；mutation_frequency为突变频率；distinct_variant_count为独立变异数；mutation_record_count为原始记录数；dominant_variant_classification及其计数为主要功能分类；variant_classification_profile和variant_type_profile为完整构成。
使用位置：补充材料中的全基因突变概况，以及正文Top mutated genes图表的源数据。
注意：mutation_frequency以421个样本为分母；如果正文最终分析人群不是421例，必须调整分母或在表注中说明。


Supplementary Table 2.csv
名称建议：Gene cell-type-specificity scores
记录数：256,207
用途：记录每个基因在不同细胞类型中的表达特异性，一行代表一个gene–Celltype组合。
主要字段：gene和Celltype为基因及细胞类型；Gene_CTS为处理后的特异性评分；Gene_CTS_original为原始评分；DetectionFraction为相应细胞类型中检测到该基因的细胞比例。
使用位置：单细胞转录组细胞类型特异性分析、CTS计算及调控模块评分的底层数据。


Supplementary Table 3.csv
名称建议：Peak accessibility across cell types
记录数：443,513
用途：记录ATAC-seq调控峰在不同细胞类型中的可及性和峰特异性，一行代表一个peak–Celltype组合。
主要字段：peak_id为峰编号或基因组坐标；Celltype为细胞类型；PeakAccessibility为峰可及性；Peak_CTS为峰的细胞类型特异性评分。
使用位置：单细胞ATAC分析、细胞类型特异性调控峰识别和Peak_CTS计算的源数据。


Supplementary Table 4.csv
名称建议：Mutation-to-regulatory-peak mapping
记录数：802,273
用途：记录体细胞变异与调控峰的空间重叠或关联关系，一行代表一个variant–peak映射。
主要字段：variant_id和peak_id为变异及峰编号；chr/start/end/ref/alt为坐标及等位基因；n_samples为涉及样本数；impact_class、func、exonic_func和gene_anno为功能注释；VariantMutationRelevance为变异与调控峰的相关性或贡献评分。
使用位置：将突变层连接到染色质可及性层的底层数据，也是MutationPeakScore和调控网络突变权重的依据。
注意：表中保留了流程连接后产生的.x和.y字段。正式公开前建议说明两组字段来源，或重命名为variant_*和matched_*，避免读者混淆。


Supplementary Table 5.csv
名称建议：Regulatory-peak mutation scores
记录数：62,076
用途：汇总每个调控峰的突变负荷和标准化突变相关评分，一行代表一个peak。
主要字段：peak_id为峰编号；n_variants为关联变异数；n_mutated_samples为突变样本数；mutation_peak_raw及capped版本为峰突变分数；MutationPeakScore为标准化评分；MutationRelevanceRaw、MutationRelevanceCapped和MutationRelevance为突变相关性分数。
使用位置：峰层面的突变评分结果，以及后续TF–peak–gene网络中MGRS_edge计算的输入。


Supplementary Table 6.csv
名称建议：Mutation-weighted TF–peak–gene regulatory network
记录数：539,403
用途：提供完整的突变加权TF–peak–gene调控边，一行代表一个特定细胞类型中的TF–peak–target gene关系。
主要字段：TF、Celltype、gene、peak_id和motif_id为调控边组成；distance为peak与靶基因距离；PeakGeneCor、PeakGeneP、PeakGeneFDR和PeakGeneStrength为peak–gene关联统计量；Target_CTS、TF_CTS、PeakAccessibility及Peak_CTS为特异性信息；CellEdgeScore和EdgeQuantile为调控边得分；n_variants、n_mutated_samples、MutationPeakScore和MutationRelevance为突变信息；MGRS_edge为该边对突变相关调控评分的贡献。
使用位置：网络图、E2F3调控模块提取、靶基因列表生成，以及MGRS/MRS计算的完整源数据。
注意：该文件约235 MB，属于机器可读Source Data，不适合直接阅读；投稿平台若限制单文件大小，可压缩为.csv.gz后上传。


Supplementary Table 7.csv
名称建议：Adjusted single-gene Cox regression results
记录数：19
用途：评价候选基因表达与总体生存的关系，并校正年龄、性别、肿瘤分期和吸烟状态。
主要字段：cohort和model为队列及模型；gene和term为基因及模型项；HR、lower95、upper95和pvalue为校正后结果；n和events为样本数及事件数；PH_gene_p和PH_global_p为比例风险假设检验；FDR为多重检验校正值。
使用位置：单基因Cox森林图和六基因候选筛选的统计依据。
注意：这是对单个基因分别建立并加入临床协变量的校正Cox模型，不等同于把全部基因同时放入一个多因素Cox模型。


Supplementary Table 8.csv
名称建议：Sample-level survival data for prognostic validation
记录数：1,497
用途：提供TCGA-LUAD及正文所用GEO队列的样本级生存信息，一行代表一个样本。
主要字段：cohort和platform_id为队列及平台；sample_name和patient_id为样本及患者编号；requested_endpoint和analysis_endpoint为计划及实际终点；survival_time_months和survival_status为统一后的生存数据；time、time_unit和status为模型字段；age_years为可获得时的年龄。
使用位置：Figure 6生存验证的样本级数据说明，以及各队列样本量和生存终点的核查。
注意：该表不含六基因risk_score和risk_group；风险评分及分组结果应以Table 4和相应风险评分源文件为准。


三、缩写说明
------------

CTS：cell-type specificity，细胞类型特异性。
MGRS：mutation-associated gene regulatory score，突变相关基因调控评分。
MRS/ctMRS：本文定义的突变调控模块评分；正文、图注和Methods中应统一最终名称。
TF：transcription factor，转录因子。
OS：overall survival，总体生存。
HR：hazard ratio，风险比。
FDR：false discovery rate，多重检验校正后的错误发现率。
PH：proportional hazards，比例风险假设。


四、推荐使用顺序
----------------

正文优先引用：Table 1 → Table 2 → Table 3 → Table 4。
突变和单细胞底层数据：Supplementary Tables 1–5。
完整调控网络：Supplementary Table 6。
预后统计与样本级验证数据：Supplementary Tables 7–8。

生成说明：本文件根据当前目录内CSV文件的实际字段和记录数编写。若后续重新筛选样本、修改评分公式或更换生存终点，应同步更新相应表格和本说明。
