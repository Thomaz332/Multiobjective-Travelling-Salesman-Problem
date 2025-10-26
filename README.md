# 🧭 Multiobjective Travelling Salesman Problem (NSGA-II)

Este projeto implementa uma versão **multiobjetivo** do clássico **Problema do Caixeiro-Viajante (TSP)**, utilizando o algoritmo genético **NSGA-II (Non-dominated Sorting Genetic Algorithm II)**.  
O objetivo é otimizar simultaneamente **três funções de custo**:  
- 🛣️ **Distância total percorrida**  
- ⏱️ **Tempo total de viagem**  
- 💰 **Custo total de pedágios**

---

## 🚀 Características principais

- Implementação **completa do NSGA-II** (ordenação não dominada e crowding distance).  
- **50 cidades** com coordenadas fixas em um plano 2D.  
- **Distâncias euclidianas** calculadas entre as cidades.  
- **Tempos de viagem** proporcionais às distâncias, com pequeno ruído aleatório para diversidade.  
- **Pedágios** aleatórios associados a cada cidade.  
- Geração da **Fronteira de Pareto** entre os objetivos.
