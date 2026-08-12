# epiPortrait 发表级 Benchmark 蓝图(综合版 v2)
## 基于 2026-08-11 两份蓝图核实后的合并与修正

**日期**:2026-08-11
**版本**:v2(合并 `benchmark_dataset_blueprint` + `tumor_benchmark_dataset_strategy`,修正数据引用与设计缺口)

---

# 1. 定位

epiPortrait = **Replicate-Aware Epigenomic Domain Profiling**

四 Claim:

- **Claim 1**:epiPortrait 将 domain 的 signal magnitude 与 native breadth geometry 分开刻画(Intensity × Breadth 分解,SignalDispersion 保留内部 architecture)。
- **Claim 2**:replicate-level call → support aggregation 减少 single-replicate-driven calls。
- **Claim 3**:检测 condition-specific domain remodeling(而非静态 Super-domain)。
- **Claim 4**:泛化到 active/repressive marks、ChIP-seq/CUT&Tag,不同 mark 用不同生物学解释。

---

# 2. 核心数据集(已核实,全部真实可获取)

## 主数据集(肿瘤优先,支撑主论文)

| 优先级 | 数据集 | 核实状态 | 设计 | 主 marks | 对应 Claim |
|---|---|---|---|---|---|
| Flagship | **GSE259262** | ✅ GBM,10 GIC vs 10 syngeneic iNSC | H3K4me3/H3K27ac/H3K27me3/H3K36me3,100 samples | 1+4 |
| 主验证 | **GSE136128** | ✅ C4-2 ENZ-R,3v3 | H3K27ac+AR+FOXA1 | 2+3 |
| 机制验证 | **GSE126004** | ✅ T47D IL6 ± FOXA1 | H3K27ac+ER/FOXA1/STAT3 | 3 |

## 补充数据集

| 优先级 | 数据集 | 核实状态 | 用途 |
|---|---|---|---|
| Supplement | **GSE78206** | ✅ MM1S ± KDM5-C70,H3K4me3 | Intensity ≠ Breadth 分离 |
| Reference | **GSE151556** | ✅ broad-H3K4me3 癌基因参考 | known-positive validation |
| Reference | **GSE145938** | ✅ MM H3K27ac SE | ROSE SE 分解 |

## 补充(非肿瘤/外部验证,来自第一份蓝图)

| 数据集 | 核实状态 | 用途 |
|---|---|---|
| **Synthetic pipeline**(自建) | 待建 | 真值 accuracy/robustness |
| ENCODE **K562**(ENCSR000AKP/AKU) | ⚠️ 编号需在 Portal 人工确认 | 历史 benchmark |
| **GSE212265** | ✅ EVT TSC/D3/D8 H3K27ac | 动态 remodeling 外部验证 |
| **GSE201262** | ✅ MCF7 WT/MCM2-2A CUT&Tag | repressive/CUT&Tag 泛化 |
| ~~GSE78913~~ | ❌ **已剔除**:实际是 TENET(前列腺 FAIRE),非 MCF7 E2 45min | 替换待定 |

---

# 3. 设计原则(来自两份蓝图,保留全部合理处)

1. **先 Claim 后数据**:每数据集对应一个明确 methodological claim,不为凑数。
2. **Synthetic truth benchmark**(自建 pipeline,不用 example_se):唯一可算 sensitivity/specificity/F1 的层。
3. **避免 circular benchmark**:Breadth 不依赖 MACS2 width 自证;主证据 = replicate reproducibility + biological enrichment + 参数鲁棒性。
4. **ROSE 是 reference 非 truth**:报告 method concordance(bp Jaccard/rank corr),并做 SE 分解(ROSE SE → Intensity/Breadth/Dual)。
5. **参数冻结**:dev 数据(K562/MCF7/GBM)调参 → 冻结 → 外部验证(GSE212265/GSE201262)。允许 mark-aware preset(stitching/promoter exclusion),不允许按结果调 elbow/support rule。
6. **统一 preprocessing**:per-replicate native peak calling(不得用 merged/consensus 替代)+ 统一 normalized BigWig。
7. **肿瘤异质性**:patient-derived samples ≠ 技术重复;10 GIC 的 support 解释为 population recurrence。
8. **数据少而强**:主 3 + 补 3 + 外部 2 = 足够,不堆 20 数据集。

---

# 4. 本版修正与补强(相对两份原始蓝图)

## 4.1 剔除 GSE78913(数据引用错误)

原文档声称 "MCF7 H3K27ac, Estradiol 45min vs EtOH",实际该编号为 **TENET**(FAIRE/ChIP in prostate/breast + C42B RNA-seq),与描述不符。**已剔除**。若需 acute E2 激活验证,改用其他真实数据集(待确认)。

## 4.2 补 SignalDispersion 可验证 metric(第一份蓝图缺口)

SignalDispersion 定位:辅助描述,不升级为第三 Super axis。但必须有直接验证:

> **同 width、同 AUC、不同内部 architecture(single compact / two-peak / plateau / fragmented)下,SignalDispersion 是否可区分?**

模拟 → 报告 SignalDispersion 分布差异(effect size + CI)。这是它唯一的增量价值,必须直接测。

## 4.3 明确 Uncertain(no-call)在 metrics 中的处理(第一份蓝图缺口)

- Precision/Recall/F1:Uncertain 域**排除**在分子分母外,单独报告 no-call rate。
- 增加 certainty calibration:高置信(如 3/3 support)call 的正确率是否高于低置信。
- 所有数据集统一报告:unanimous rate / majority-supported rate / no-call rate(呼应 Claim 2)。

## 4.4 补 per-replicate native peak 上游工具链(第一份蓝图缺口)

Breadth 依赖 per-replicate native peaks。统一建议:
- **caller**:MACS2(标准;broad 场景可用 `--broad`),或 CUT&Tag 适配 caller;
- 每个 replicate 独立 call(禁 merged/consensus 替代);
- 产出 `peak_path` → epiPortrait `feature="Breadth"`。

---

# 5. Benchmark 层与 Metric 矩阵

## 5.1 Synthetic truth(自建)

- 5,000–20,000 domains,2 conditions,3 rep/condition
- 真值:Typical/Intensity-only/Breadth-only/Dual
- 两套比例:balanced(25/25/25/25)+ realistic(70/10/10/10)
- 干扰:replicate CV(10/20/30%)、boundary jitter(±5/10/20%)、dropout、one-rep outlier
- Metric:Intensity/Breadth 各算 Sens/Spec/Prec/Recall/F1;Combined:Macro-F1、Cohen's kappa、accuracy;外加鲁棒性 3 项

## 5.2 H3K4me3 Breadth(GSE259262 + K562/MCF7)

- literature concordance(MACS2 broad top5% 仅作参考,不作 truth)
- replicate/population support rate、Cohen's kappa
- Breadth-Super 富集 cell identity / lineage gene sets(odds ratio + 95%CI)
- matched RNA(如有):median expression + log2 diff

## 5.3 H3K27ac ROSE(GSE259262 + GSE145938 + K562)

- ROSE 同预处理(candidate/stitching/promoter exclusion 统一)
- bp Jaccard / region overlap / Spearman rank / ROSE SE recall
- **ROSE SE 分解**:ROSE SE → Intensity/Breadth/Dual,比较 signal/width/SignalDispersion/expression
- 外部:rE2G/CRISPR regulatory-link enrichment(只作 enrichment,不 claim 因果 target)

## 5.4 Dynamic remodeling(GSE136128 + GSE126004)

- per_group → support → transition(Typical→Intensity-Super 等)
- 报告 Intensity/NativeMaxPeakWidth/SignalDispersion 变化
- GSE136128:新获得 domain 是否伴随 AR/FOXA1 occupancy gain(orthogonal)
- GSE126004:FOXA1-dependent vs independent transitions

## 5.5 Repressive/CUT&Tag 泛化(GSE201262)

- H3K27me3/H3K9me3,mark-aware 术语(Extended-Domain/Intensity-Extreme/Dual-Extreme)
- width expansion/contraction、Intensity change、combined remodeling
- RNA 表达 + 原始研究报道区域验证

## 5.6 参数鲁棒性(对所有主数据)

- shared-domain universe 两套(conservative/default)Jaccard
- min_peak_overlap_fraction 0.25/0.50/0.75
- peak caller 两套(MACS2 standard/broad)
- leave-one-replicate-out(validation script,非 API)
- runtime:1k/10k/50k/100k domains × 4/6/12 samples

---

# 6. 统一输出(所有数据集)

```text
replicate call concordance
support fraction
class composition
Intensity distribution
NativeMaxPeakWidth
SignalDispersion
Combined class
transition
```

比较类数据集追加:bp Jaccard、Cohen's kappa、Spearman rho、odds ratio、95% CI。

统计报告原则:effect size + 95% CI + distribution,不只 p-value。

---

# 7. 推荐 Figure

- Fig1:framework
- Fig2:Synthetic truth + robustness
- Fig3:GBM 多 mark(GSE259262)
- Fig4:前列腺耐药(GSE136128)+ AR/FOXA1 验证
- Fig5:乳腺癌 IL6/FOXA1(GSE126004)
- Fig6:跨癌种 domain phenotype 对比
- FigS:H3K4me3 perturbation(GSE78206)、broad-H3K4me3 reference(GSE151556)、MM SE 分解(GSE145938)、repressive(CUT&Tag GSE201262)、参数鲁棒性、runtime

---

# 8. 执行优先级

## 第一批(立刻)
1. Synthetic truth pipeline(R 实现)
2. GSE259262(flagship 多 mark)
3. GSE136128(耐药 remodeling)

## 第二批
4. GSE126004(机制验证)
5. GSE78206(Intensity/Breadth 分离)

## 第三批(reference validation)
6. GSE151556 / GSE145938

## 外部验证(参数冻结后)
7. GSE212265 / GSE201262

---

# 9. 一句话结论

> **以 Synthetic truth → GBM 多标记(GSE259262)→ 前列腺耐药(GSE136128)→ 乳腺癌信号重塑(GSE126004)为证据链主线,补以 KDM5/MM 的 Intensity/Breadth 分离(GSE78206)、已知 broad/SE 参考(GSE151556/GSE145938),再用 trophoblast 与 repressive CUT&Tag 做冻结参数后的外部验证。数据少而强,每集解决一个明确 Claim,经得起 Reviewer 推敲。**
