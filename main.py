from AGFuncs import *
import numpy as np
from typing import Final
import matplotlib.pyplot as plt
import functools
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

# --- PARÂMETROS DO ALGORITMO ---
TAM_CROMO: Final[int] = len(CIDADES)
TAM_POP: Final[int] = TAM_CROMO * 10
NUM_GERACOES: Final[int] = 100
TAXA_MUTACAO: Final[float] = 0.3
TAXA_CRUZAMENTO: Final[float] = 0.8
contador_estagnacao = 0
GERACOES_PARA_PARAR = 50 
qtd_geracoes = 0

if __name__ == '__main__':

    populacao = np.zeros((TAM_POP, TAM_CROMO), dtype=int)
    nota_populacao = np.zeros((TAM_POP, 6))
    iniciaPopulacao(populacao)
    print("--- POPULAÇÃO INICIAL GERADA ---")
    
    calculoNotas(populacao, nota_populacao)
    
    historico_plot = []

    with open('historico_fronteiras.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Geracao', 'ID_Original', 'Cromossomo', 'Distancia_KM', 'Tempo_H', 'Pedagio_R$'])

        for geracao in range(NUM_GERACOES):
            qtd_geracoes = geracao
            
            ranks, fronts = calculoFronteDePareto(nota_populacao, TAM_POP)
            crowding_distances = calculate_crowding_metrics(nota_populacao, fronts)
            if fronts[0] == 100:
                break
            if (geracao + 1) % 10 == 0 or geracao == 0:
                 print(f"\rGeração {geracao+1}/{NUM_GERACOES} | Membros na Fronteira: {len(fronts[0])}", end="")
                 print_fronteira_geracao(geracao + 1, populacao, nota_populacao, fronts)

            melhores_indices_geracao = fronts[0]
            for idx in melhores_indices_geracao:
                rota_str = " - ".join([str(int(c)) for c in populacao[idx]])
                dist = nota_populacao[idx, 2]
                tempo = nota_populacao[idx, 3]
                custo = nota_populacao[idx, 4]
                csv_writer.writerow([geracao + 1, idx, rota_str, f"{dist:.2f}", f"{tempo:.2f}", f"{custo:.2f}"])
            historico_plot.append(nota_populacao[melhores_indices_geracao, 2:5])

            populacao_filhos = np.zeros_like(populacao)
            for i in range(TAM_POP):
                idx_pai1 = selecao_NSGA2(ranks, crowding_distances)
                idx_pai2 = selecao_NSGA2(ranks, crowding_distances)
                pai1 = populacao[idx_pai1]
                pai2 = populacao[idx_pai2]
                
                if random.random() < TAXA_CRUZAMENTO:
                    filho = crossoverOX(pai1, pai2, TAM_CROMO)
                else:
                    filho = pai1.copy()
                filho = mutacao_swap(filho, TAXA_MUTACAO)
                populacao_filhos[i] = filho
            
            populacao_combinada = np.vstack([populacao, populacao_filhos])
            
            nota_populacao_combinada = np.zeros((TAM_POP * 2, 6))
            calculoNotas(populacao_combinada, nota_populacao_combinada)

            ranks_combinados, fronts_combinados = calculoFronteDePareto(nota_populacao_combinada, TAM_POP * 2)
            crowding_combinados = calculate_crowding_metrics(nota_populacao_combinada, fronts_combinados)
            if len(fronts[0]) == TAM_POP:
                contador_estagnacao += 1
            else:
                contador_estagnacao = 0

            if contador_estagnacao >= GERACOES_PARA_PARAR:
                print(f"\nParando na Geração {geracao+1} devido à estagnação da fronteira por {GERACOES_PARA_PARAR} gerações.")
                break
            proxima_populacao = np.zeros_like(populacao)
            proxima_nota_populacao = np.zeros_like(nota_populacao)
            
            prox_idx_livre = 0
            front_num = 0
            
            while prox_idx_livre < TAM_POP:
                indices_da_fronteira_atual = fronts_combinados[front_num]
                
                if prox_idx_livre + len(indices_da_fronteira_atual) > TAM_POP:
                    break 
                
                for idx in indices_da_fronteira_atual:
                    proxima_populacao[prox_idx_livre] = populacao_combinada[idx]
                    proxima_nota_populacao[prox_idx_livre] = nota_populacao_combinada[idx]
                    prox_idx_livre += 1
                front_num += 1
            
            if prox_idx_livre < TAM_POP:
                ultima_fronteira = fronts_combinados[front_num]
                
                crowding_da_fronteira = crowding_combinados[ultima_fronteira]
                
                indices_ordenados = np.argsort(-crowding_da_fronteira)
                
                for i in range(TAM_POP - prox_idx_livre):
                    idx_original_na_comb = ultima_fronteira[indices_ordenados[i]]
                    proxima_populacao[prox_idx_livre + i] = populacao_combinada[idx_original_na_comb]
                    proxima_nota_populacao[prox_idx_livre + i] = nota_populacao_combinada[idx_original_na_comb]

            populacao = proxima_populacao
            nota_populacao = proxima_nota_populacao

        print("\n--- FRONTEIRA DE PARETO DA GERAÇÃO FINAL ---")

        final_ranks, final_fronts = calculoFronteDePareto(nota_populacao, TAM_POP)
        melhores_indices_finais = final_fronts[0]

        headers = ["ID na Pop. Final", "Rota (Cromossomo)", "Distância (KM)", "Tempo (H)", "Pedágio (R$)"]
        table_data = []

        for idx in melhores_indices_finais:
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

    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(111, projection='3d')
    cores = cm.viridis(np.linspace(0, 1, len(historico_plot)))
    for i, fronteira in enumerate(historico_plot):
        if len(fronteira) > 0:
            ax.scatter(fronteira[:, 0], fronteira[:, 1], fronteira[:, 2], color=cores[i], alpha=0.4, s=15)
    fronteira_inicial = historico_plot[0]
    fronteira_final = historico_plot[-1]
    if len(fronteira_inicial) > 0:
        ax.scatter(fronteira_inicial[:, 0], fronteira_inicial[:, 1], fronteira_inicial[:, 2], 
                   color='cyan', s=150, marker='o', label='Geração 1', depthshade=False, zorder=10)
    if len(fronteira_final) > 0:
        ax.scatter(fronteira_final[:, 0], fronteira_final[:, 1], fronteira_final[:, 2], 
                   color='red', s=150, marker='X', label=f'Geração Final ({qtd_geracoes})', depthshade=False, zorder=10)
    ax.set_title('Convergência da Fronteira de Pareto 3D (NSGA-II)', fontsize=16)
    ax.set_xlabel('Distância Total (KM)', fontsize=12)
    ax.set_ylabel('Tempo Total (H)', fontsize=12)
    ax.set_zlabel('Custo Pedágio (R$)', fontsize=12)
    ax.legend()
    ax.view_init(elev=20, azim=45)
    plt.savefig('convergencia_pareto_todas_geracoes.png')
    plt.figure(figsize=(18, 12))

    # Gráfico 1: Distância vs Tempo
    plt.subplot(2, 2, 1)
    for i, fronteira in enumerate(historico_plot):
        if len(fronteira) > 0:
            plt.scatter(fronteira[:, 0], fronteira[:, 1], color=cores[i], alpha=0.4, s=15)
    if len(fronteira_inicial) > 0:
        plt.scatter(fronteira_inicial[:, 0], fronteira_inicial[:, 1], 
                color='cyan', s=150, marker='o', label='Geração 1', zorder=10)
    if len(fronteira_final) > 0:
        plt.scatter(fronteira_final[:, 0], fronteira_final[:, 1], 
                color='red', s=150, marker='X', label=f'Geração Final ({qtd_geracoes})', zorder=10)
    plt.xlabel('Distância Total (KM)')
    plt.ylabel('Tempo Total (H)')
    plt.title('Distância vs Tempo')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Gráfico 2: Distância vs Pedágio
    plt.subplot(2, 2, 2)
    for i, fronteira in enumerate(historico_plot):
        if len(fronteira) > 0:
            plt.scatter(fronteira[:, 0], fronteira[:, 2], color=cores[i], alpha=0.4, s=15)
    if len(fronteira_inicial) > 0:
        plt.scatter(fronteira_inicial[:, 0], fronteira_inicial[:, 2], 
                color='cyan', s=150, marker='o', label='Geração 1', zorder=10)
    if len(fronteira_final) > 0:
        plt.scatter(fronteira_final[:, 0], fronteira_final[:, 2], 
                color='red', s=150, marker='X', label=f'Geração Final ({qtd_geracoes})', zorder=10)
    plt.xlabel('Distância Total (KM)')
    plt.ylabel('Custo Pedágio (R$)')
    plt.title('Distância vs Pedágio')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Gráfico 3: Tempo vs Pedágio
    plt.subplot(2, 2, 3)
    for i, fronteira in enumerate(historico_plot):
        if len(fronteira) > 0:
            plt.scatter(fronteira[:, 1], fronteira[:, 2], color=cores[i], alpha=0.4, s=15)
    if len(fronteira_inicial) > 0:
        plt.scatter(fronteira_inicial[:, 1], fronteira_inicial[:, 2], 
                color='cyan', s=150, marker='o', label='Geração 1', zorder=10)
    if len(fronteira_final) > 0:
        plt.scatter(fronteira_final[:, 1], fronteira_final[:, 2], 
                color='red', s=150, marker='X', label=f'Geração Final ({qtd_geracoes})', zorder=10)
    plt.xlabel('Tempo Total (H)')
    plt.ylabel('Custo Pedágio (R$)')
    plt.title('Tempo vs Pedágio')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Gráfico 4: Evolução do tamanho da fronteira de Pareto
    plt.subplot(2, 2, 4)
    tamanhos_fronteiras = [len(fronteira) for fronteira in historico_plot]
    geracoes = range(1, len(tamanhos_fronteiras) + 1)
    plt.plot(geracoes, tamanhos_fronteiras, 'b-', linewidth=2, marker='o', markersize=4)
    plt.xlabel('Geração')
    plt.ylabel('Tamanho da Fronteira')
    plt.title('Evolução do Tamanho da Fronteira de Pareto')
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('analise_pareto_2d.png')

    print("\n--- ANÁLISE ESTATÍSTICA DA FRONTEIRA FINAL ---")
    if len(fronteira_final) > 0:
        print(f"Total de soluções na fronteira: {len(fronteira_final)}")
        print(f"Distância: Min={fronteira_final[:, 0].min():.2f}KM, Max={fronteira_final[:, 0].max():.2f}KM")
        print(f"Tempo: Min={fronteira_final[:, 1].min():.2f}H, Max={fronteira_final[:, 1].max():.2f}H")
        print(f"Pedágio: Min={fronteira_final[:, 2].min():.2f}R$, Max={fronteira_final[:, 2].max():.2f}R$")
    plt.show()