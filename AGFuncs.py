import sys
import random
import numpy as np
from typing import Final
from tabulate import tabulate

CIDADES: Final[dict] = {
    0: (5, 5), 1: (10, 12), 2: (15, 8), 3: (8, 15), 4: (90, 5),
    5: (85, 80), 6: (5, 80), 7: (45, 50), 8: (50, 45), 9: (20, 95),
    10: (25, 25), 11: (30, 35), 12: (40, 40), 13: (55, 60), 14: (65, 70),
    15: (75, 25), 16: (80, 30), 17: (35, 75), 18: (60, 20), 19: (70, 85),
    20: (15, 45), 21: (25, 55), 22: (45, 65), 23: (65, 45), 24: (85, 15),
    25: (10, 65), 26: (30, 85), 27: (50, 25), 28: (70, 65), 29: (90, 35),
    30: (20, 20), 31: (40, 30), 32: (60, 50), 33: (80, 70), 34: (95, 90),
    35: (12, 35), 36: (32, 45), 37: (52, 55), 38: (72, 75), 39: (92, 85),
    40: (18, 75), 41: (38, 65), 42: (58, 35), 43: (78, 55), 44: (98, 25),
    45: (22, 85), 46: (42, 15), 47: (62, 85), 48: (82, 5), 49: (95, 95)
}   

def iniciaPopulacao(populacao):
    # Cria o conjunto de população com valores aleatorios
    for i in range(len(populacao)):
        populacao[i] = np.random.permutation(len(CIDADES))

# Calculo feito da distancia entre 2 cidades, cidadem1 e 2, depois eleva ambos ao quadrado pra retirar qual

def distancia(cidade1, cidade2):
    x1, y1 = CIDADES[(cidade1)]
    x2, y2 = CIDADES[(cidade2)]
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

TEMPO = np.zeros((len(CIDADES), len(CIDADES)))
for i in range(len(CIDADES)):
    for j in range(len(CIDADES)):
        if i == j:
            TEMPO[i][j] = None
        else:
            dist = distancia(i, j)
            ruido = np.random.uniform(0.06, 0.7)
            TEMPO[i][j] = dist * ruido

PRECO_PEDAGIO = np.zeros((len(CIDADES), len(CIDADES)))

for i in range(len(CIDADES)):
    for j in range(len(CIDADES)):
        if i == j:
            PRECO_PEDAGIO[i][j] = None
        else:
            dist = distancia(i, j)
            fator_aleatorio = np.random.uniform(0.06, 0.4)
            PRECO_PEDAGIO[i][j] = max(5, min(150, dist * 0.5 * fator_aleatorio))

# Distancia total do percurso de cada cromossomo

def distanciaTotal(cromossomo):
    total = 0
    for i in range(len(cromossomo) - 1):
        total += distancia(cromossomo[i], cromossomo[i+1])
    # Volta para a cidade inicial
    total += distancia(cromossomo[-1], cromossomo[0])
    return total


def tempoTotal(cromossomo):
    total = 0
    for i in range(len(cromossomo) - 1):
        total += TEMPO[int(cromossomo[i])][int(cromossomo[i+1])]
    # Volta à cidade inicial
    total += TEMPO[int(cromossomo[-1])][int(cromossomo[0])]
    return total


def custoPedagio(cromossomo):
    total = 0
    for i in range(len(cromossomo) - 1):
        total += PRECO_PEDAGIO[int(cromossomo[i])][int(cromossomo[i+1])]
    
    total += PRECO_PEDAGIO[int(cromossomo[-1])][int(cromossomo[0])]
    return total

#       --NOTA POPULACAO--
# [i][0] -> Indice original na matriz 
# [i][1] -> Nota 
# [i][2] -> Distancia total 
# [i][3] -> Tempo total
# [i][4] -> preco pedágio
# [i][5] -> percentual de seleção

# Fitness de cada cromossomo, por enquanto apenas distancia

def calculoNotas(populacao, nota_populacao):
    for i in range(len(populacao)):
        nota_populacao[i, 0] = i
        nota_populacao[i, 2] = distanciaTotal(populacao[i])
        nota_populacao[i, 3] = tempoTotal(populacao[i])
        nota_populacao[i, 4] = custoPedagio(populacao[i])


def domina(indice_a, indice_b, nota_populacao):
    objetivos_a = nota_populacao[indice_a, 2:5]
    objetivos_b = nota_populacao[indice_b, 2:5]
    return np.all(objetivos_a <= objetivos_b) and np.any(objetivos_a < objetivos_b)


def calculoFronteDePareto(nota_populacao, TAM_POP):
    domination_counts = np.zeros(TAM_POP, dtype=int)
    dominated_solutions = [[] for _ in range(TAM_POP)]
    ranks = np.zeros(TAM_POP, dtype=int)
    fronts = [[]]

    for i in range(TAM_POP):
        for j in range(i + 1, TAM_POP):
            if domina(i, j, nota_populacao):
                dominated_solutions[i].append(j)
                domination_counts[j] += 1
            elif domina(j, i, nota_populacao):
                dominated_solutions[j].append(i)
                domination_counts[i] += 1
    
    for i in range(TAM_POP):
        if domination_counts[i] == 0:
            ranks[i] = 0
            fronts[0].append(i)
            
    rank_atual = 0
    while rank_atual < len(fronts):
        next_front = []
        for i in fronts[rank_atual]:
            for j in dominated_solutions[i]:
                domination_counts[j] -= 1
                if domination_counts[j] == 0:
                    ranks[j] = rank_atual + 1
                    next_front.append(j)
        rank_atual += 1
        if len(next_front) > 0:
            fronts.append(next_front)
            
    return ranks, fronts


def calculate_crowding_metrics(nota_populacao, fronts):
    num_individuals = len(nota_populacao)
    crowding_metrics = np.zeros(num_individuals)
    
    objetivos = nota_populacao[:, 2:5]
    num_objectives = objetivos.shape[1]

    for front in fronts:
        if not front: continue
        
        front_objetivos = objetivos[front, :]
        
        for m in range(num_objectives):
            sorted_indices = np.argsort(front_objetivos[:, m])
            sorted_front = np.array(front)[sorted_indices]
            
            crowding_metrics[sorted_front[0]] = np.inf
            crowding_metrics[sorted_front[-1]] = np.inf
            
            if len(sorted_front) > 2:
                min_obj = front_objetivos[sorted_indices[0], m]
                max_obj = front_objetivos[sorted_indices[-1], m]
                range_obj = max_obj - min_obj
                if range_obj == 0: continue

                for i in range(1, len(sorted_front) - 1):
                    dist = objetivos[sorted_front[i+1], m] - objetivos[sorted_front[i-1], m]
                    crowding_metrics[sorted_front[i]] += dist / range_obj
                    
    return crowding_metrics


def selecao_NSGA2(ranks, crowding_metrics, k=2):
    pop_size = len(ranks)
    
    # Sorteia k competidores (k=2 para torneio binário)
    candidatos_indices = np.random.choice(pop_size, k, replace=False)
    
    idx1 = candidatos_indices[0]
    idx2 = candidatos_indices[1]
    
    if ranks[idx1] < ranks[idx2]:
        return idx1
    elif ranks[idx2] < ranks[idx1]:
        return idx2
    
    if crowding_metrics[idx1] > crowding_metrics[idx2]:
        return idx1
    else:
        return idx2


def mutacao_swap(cromossomo, TAXA_MUTACAO):
    if random.random() < TAXA_MUTACAO:
        idx1, idx2 = np.random.choice(len(cromossomo), 2, replace=False)
        cromossomo[idx1], cromossomo[idx2] = cromossomo[idx2], cromossomo[idx1]
    return cromossomo


def crossoverOX(pai1, pai2, TAM_CROMO):
    filho = np.full(TAM_CROMO, -1)  

    start, end = sorted(np.random.randint(0, len(CIDADES), 2))
    
    filho[start:end+1] = pai1[start:end+1]
    
    pos_filho = (end+1) % len(CIDADES)
    for c in pai2:
        if c not in filho:
            filho[pos_filho] = c
            pos_filho = (pos_filho + 1) % len(CIDADES)
    
    return filho

def printPopulacao(populacao, notaPopulacao):
    headers = ["ID", "CROMOSSOMO", "DISTÂNCIA", "TEMPO", "PEDÁGIO"]
    table_data = []

    for i in range(len(populacao)):
        rota_str = " - ".join([str(int(cidade)) for cidade in populacao[i]])
        distancia = notaPopulacao[i][2]
        tempo = notaPopulacao[i][3]
        custo = notaPopulacao[i][4]
        table_data.append([i, rota_str, f"{distancia:.2f}KM", f"{tempo:.2f}H", f"{custo:.2f}R$"])

    print("\nPOPULAÇÃO ATUAL:")
    print(tabulate(table_data, headers=headers, tablefmt="grid"))

def print_fronteira_geracao(geracao, populacao, nota_populacao, fronts):
    print(f"\n\n--- FRONTEIRA DE PARETO | GERAÇÃO {geracao} ---")

    melhores_indices = fronts[0]
    
    headers = ["ID na Pop.", "Rota (Cromossomo)", "Distância (KM)", "Tempo (H)", "Pedágio (R$)"]
    table_data = []

    for idx in melhores_indices:
        cromossomo = populacao[idx]
        rota_str = " -> ".join(map(str, cromossomo[:10])) + f"... -> {cromossomo[0]}"
        
        dist = nota_populacao[idx, 2]
        tempo = nota_populacao[idx, 3]
        custo = nota_populacao[idx, 4]
        
        table_data.append([
            idx,
            rota_str,
            f"{dist:.2f}",
            f"{tempo:.2f}",
            f"{custo:.2f}"
        ])

    print(tabulate(table_data, headers=headers, tablefmt="grid"))
    print()