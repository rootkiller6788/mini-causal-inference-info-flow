# Mini 因果推断与信息流

**从零开始、零依赖的 C 语言实现**，覆盖因果推断理论与信息论因果分析。每个模块将 Pearl 结构因果模型、do-演算、反事实推理、Granger 因果以及 Schreiber 传递熵等理论方程翻译为可运行的 C 代码。

## 子模块

| 子模块 | 主题 | 核心课程 |
|--------|------|----------|
| [mini-causal-discovery-pc-fci](mini-causal-discovery-pc-fci/) | PC/FCI 算法、条件独立性检验、基于约束的因果结构发现、Bootstrap 稳定性评估 | Stanford STATS 366, CMU 10-708 |
| [mini-causal-effect-estimation-dag](mini-causal-effect-estimation-dag/) | SCM、DAG 协变量调整集、后门准则、前门准则、因果效应识别 | Harvard EPI 289, MIT 6.438 |
| [mini-counterfactual-reasoning](mini-counterfactual-reasoning/) | 溯因-行动-预测三步法、反事实界限、潜在结果、中介分析 | Harvard EPI 289, Stanford MS&E 338 |
| [mini-do-calculus-intervention](mini-do-calculus-intervention/) | do-演算三规则、图手术、干预、可识别性、因果效应估计 | Stanford CS 229M, Harvard EPI 289 |
| [mini-granger-causality-time-series](mini-granger-causality-time-series/) | Granger F 检验、谱/非线性/条件 Granger、预测误差方差分解、VAR 因果 | MIT 14.384, LSE EC484 |
| [mini-info-flow-in-networks](mini-info-flow-in-networks/) | 传递熵、Granger 因果、互信息在网络拓扑上的量化、因果网络重建 | MIT 6.441, MIT 6.207 |
| [mini-pearl-structural-causal-models](mini-pearl-structural-causal-models/) | SCM 形式化、可迁移性、选择偏差修正、g-公式、反事实公平性 | Stanford STATS 362, Harvard EPI 289 |
| [mini-transfer-entropy-schreiber](mini-transfer-entropy-schreiber/) | Schreiber 传递熵、信息论因果度量、嵌入、统计显著性检验 | MIT 6.441, MPI-PKS |

## 设计哲学

- **零外部依赖** — 纯 C（C99/C11），仅依赖 `libc` 和 `libm`
- **模块自包含** — 每个目录均有独立的 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码映射** — 每个模块直接对应经典文献（Pearl 2009、Spirtes et al. 2000、Schreiber 2000、Granger 1969）
- **实用演示** — 从 CSV 数据中学习因果图、do-演算可识别性检查器、Granger 因果显著性检验、传递熵估计器等

## 构建

每个模块独立构建。进入模块目录并运行：

```bash
cd mini-causal-discovery-pc-fci
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-causal-inference-info-flow/
├── mini-causal-discovery-pc-fci/       # PC 与 FCI 基于约束的因果发现
├── mini-causal-effect-estimation-dag/  # 基于 DAG 的因果效应估计与协变量调整
├── mini-counterfactual-reasoning/      # Pearl 三步反事实推理
├── mini-do-calculus-intervention/      # do-演算、干预与图手术
├── mini-granger-causality-time-series/ # 时间序列 Granger 因果分析
├── mini-info-flow-in-networks/         # 网络上的信息流与因果分析
├── mini-pearl-structural-causal-models/# Pearl 结构因果模型与高级主题
└── mini-transfer-entropy-schreiber/    # Schreiber 传递熵与信息传递
```

## 许可证

MIT
