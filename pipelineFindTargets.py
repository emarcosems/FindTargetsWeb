from __future__ import print_function

import traceback
from datetime import datetime
from bioservices import KEGG
from bioservices import UniProt
from bioservices import NCBIblast
from bs4 import BeautifulSoup
from zipfile import ZipFile, ZIP_DEFLATED
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
from libsbml import *
from scipy.sparse import csr_matrix
from scipy.linalg import pinv
from decimal import Decimal, getcontext

from cobra.flux_analysis import single_reaction_deletion




import cobra.io.sbml3

import retrieve
import math
import xlwt
import requests
import os
import re
import sys
import cobra.manipulation
import pandas as pd
import time
import threading
import smtplib
from requests import Response
import xml.etree.ElementTree as ET
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import csv


from scipy.linalg import pinv


# CLASSE CRIADA COM O OBJETIVO DE MANTER A CHAMADA ASSINCRONA DO METODO DA MODELAGEM
class MyThread(threading.Thread):

    name = None
    organism = None
    email = None
    model = None
    method = None

    def __init__(self, name, organism, email, modelParam, method):
        threading.Thread.__init__(self)
        self.model = modelParam
        self.organism = organism
        self.name = name
        self.email = email
        self.method = method

    def run(self):
        print("A THREAD COMECOU!!!!!!!!")
        FindTargets().mainMethod(self.model, self.organism, self.name, self.email, self.method)

class FindTargets:

    # ADD ALL ORGANISMS MAPPED FROM KEGG TO ATTRIBUTE IN COMBOBOX FORM
    # USED TO FULL OUR COMBOBOX FROM INDEX APP
    def list_organism(self):
        dir_path = os.path.dirname(os.path.realpath(__file__)) # identifica o local real onde esse arquivo esta
        with open(dir_path+"/organismos.txt", "r") as infile:
            data = infile.read()
        my_list = data.splitlines()
        my_list_aux = []
        my_list_aux.append(("", " --- SELECT --- "))

        for itens in my_list:
            itens_splt = itens.split(",")
            tuple_add = (itens_splt[0], itens_splt[1])
            my_list_aux.append(tuple_add)

        return tuple(my_list_aux)

    # METODO PARA MONTAR O OBJETO COBRA PARA TRAFEGO NA APLICACAO!

    def readModel(self, sbmlfile):
        return cobra.io.sbml3.read_sbml_model(sbmlfile)

    # MODEL VALIDATION METHOD TO CONTINUE OUR EXECUTION

    def validateModel(self, model):
        return_dict_val_model = {}

        numCasasDec = 6

        initialSolution = model.optimize()
        valInitialSolutionFmt = round(float(initialSolution.objective_value), numCasasDec)
        return_dict_val_model['valInitialSolutionFmt'] = valInitialSolutionFmt

        if initialSolution.status != 'optimal':
            return_dict_val_model['message'] = "ERROR! SELECTED MODEL HAS AN ERROR. PLEASE VERIFY AND TRY AGAIN!"
            return_dict_val_model['ehParaFazer'] = False
        else:
            if valInitialSolutionFmt == 0.0 or valInitialSolutionFmt == -0.0:
                return_dict_val_model['message'] = "SELECTED MODEL DOES NOT GENERATE BIOMASS. PLEASE VERIFY AND TRY AGAIN!"
                return_dict_val_model['ehParaFazer'] = False
            else:
                return_dict_val_model['message'] = "SELECTED MODEL GENERATES BIOMASS. PRESS 'SUBMIT' TO CONTINUE!"
                return_dict_val_model['ehParaFazer'] = True

        return return_dict_val_model

    #
    # MAIN METHOD THAT REALIZE ANALYSIS OF SELECTED NETWORK SBML FILE
    #

    def mainMethod(self, model, organismParam, name, email, method):

        numCasasDec = 6
        TIMEOUT_SECONDS = 500000
        model = model

        # BEFORE, WE NEED TO CREATE FOLDERS TO INCLUDE OUR FILES AND AFTER ZIP THEM
        dataHoraAtual = datetime.now()
        dataHoraAtualFmt = dataHoraAtual.strftime('%Y-%m-%d_%H.%M.%S.%f')
        dir_path = os.path.dirname(os.path.realpath(__file__))
        directory = dir_path+"/results/"+dataHoraAtualFmt
        static_dir = dir_path.replace('/FindTargetsWEB', '/static')

        # TESTS LINES!
        #directory = dir_path+"/results/"+"testesMeriguetiCCBH"
        #directory = dir_path+"/results/"+"testesMeriguetiPAO1"
        #directory = dir_path+"/results/"+"testesMeriguetiPAO1_2017"

        dir_blasts = directory+"/blasts"

        if not os.path.exists(directory):
            os.makedirs(directory)

        if not os.path.exists(dir_blasts):
            os.makedirs(dir_blasts)

        start = time.time()
        fileLog = open(directory+"//"+"LOG_EXECUCAO_"+str(dataHoraAtualFmt)+".txt", "w")
        fileLog.write('*** INICIO DA EXECUCAO *** \n')
        fileLog.write('Modelo/Organismo/Nome/Email/Metodo = {0}, {1}, {2}, {3}, {4}'.format(model, organismParam, name, email, method))
        fileLog.write('\n\n\n')
        try:
            ####################################################
            # ITEM 1 - ANALYSIS OF P. AERUGINOSA MODEL SELECTED
            ####################################################
            fileLog.write('INICIO DO ITEM 1\n')
            self.reportModel(model, directory+"//model_data_"+str(model)+"_.xlsx")
            fileLog.write('GRAVOU ARQUIVO COM OS DADOS DE GENE/REACAO DA REDE\n')

            initialSolution = model.optimize()
            fileLog.write('GEROU FBA \n')
            fileLog.write('FBA = {0}\n'.format(initialSolution))
            valInitialSolutionFmt = round(float(initialSolution.objective_value), numCasasDec)
            
            ################################################
            # ITEM 2 - PRIORITIZATION OF EXISTING REACTIONS
            ################################################
            fileLog.write('INICIO DO ITEM 2\n')
            fva_result = cobra.flux_analysis.flux_variability_analysis(model, model.reactions[:len(model.reactions)])
            pd.DataFrame.from_dict(fva_result).round(numCasasDec).to_csv(directory+"/01-fva.csv")
            fileLog.write('GEROU FVA \n')
            listReactPerturb = []
            reacoesFVADifZero = open(directory+"//"+"02-reacoesFVADifZero.txt", "w")

            df_fva_result = pd.DataFrame(fva_result)
            df_fva_result = df_fva_result.sort_index()
            itemsFVAResult = df_fva_result.iterrows()
            itemsFBAResult = sorted(initialSolution.fluxes.items())
            fileLog.write('INICIO DA VERIFICACAO ENTRE FBA/FVA \n')
           
            
            for (keyFVA, valueFVA), (keyFBA, valueFBA) in list(zip(itemsFVAResult, itemsFBAResult)):
                minFVA = round(float(valueFVA["minimum"]), numCasasDec)
                maxFVA = round(float(valueFVA["maximum"]), numCasasDec)
                valueFBA = round(float(valueFBA), numCasasDec)

                # Filtrando reacoes que estao ativas no momento
                if maxFVA != 0:
                    # method --> ("1", "FBA+FVA"), ("2", "Only FBA")
                    if method == "1":
                        diferencaReacao = maxFVA - minFVA
                        if diferencaReacao == 0:
                            reacoesFVADifZero.write(str(keyFBA) + '\n')
                            listReactPerturb.append(keyFBA)
                    else:
                        reacoesFVADifZero.write(str(keyFBA) + '\n')
                        listReactPerturb.append(keyFBA)
                    
                    #diferencaReacao = maxFVA - minFVA
                    #reacoesFVADifZero.write("Resultado da diferenca entre max e min = " + str(diferencaReacao))
                    #reacoesFVADifZero.write('\n\n')

                if keyFBA != keyFVA:
                    raise Exception("Por favor, verifique, chave FBA/FVA nao batem na comparacao")

                else:
                    if valueFBA < minFVA:
                        #print(keyFBA)
                        #print(valueFBA)
                        #print(minFVA)
                        #print(maxFVA)
                        raise Exception("Por favor, verifique o modelo. Valor do FBA dessa reacao ser menor que o minimo descrito no FVA")

                    if valueFBA > maxFVA:
                        #print(keyFBA)
                        #print(valueFBA)
                        #print(minFVA)
                        #print(maxFVA)
                        raise Exception("Por favor, verifique o modelo. Valor do FBA dessa reacao ser maior que o maximo descrito no FVA")

            reacoesFVADifZero.close()
            fileLog.write("Lista de reacoes a analisar = " + str(len(listReactPerturb)) + '\n')
            
            
            
            fileLog.write("GENERATED FILE 01-fva.csv\n")
            fileLog.write("GENERATED FILE %s\n" % reacoesFVADifZero.name)
            
            ################################################
            # ITEM 2.1 - FLUX BALANCE ANALYSIS (FBA)
            ################################################
            
            fileLog.write('INICIO DO ITEM 2.1 - FLUX BALANCE ANALYSIS (FBA)\n')

		# 1. Definição do Sistema: Representação da Matriz Estequiométrica (S)
		# A matriz estequiométrica (S) é uma representação matemática das reações metabólicas.
		# Cada linha representa um metabólito, e cada coluna representa uma reação.
		# Os coeficientes estequiométricos indicam a participação dos metabólitos nas reações.

		# Construção explícita da matriz estequiométrica (S)
            data = []
            row_ind = []
            col_ind = []
            metabolitos = model.metabolites
            reacoes = model.reactions
            
            
            # Itera sobre cada metabólito e cada reação para construir a matriz S
            for i, metabolito in enumerate(metabolitos):
                for j, reacao in enumerate(reacoes):
                    # Obtém o coeficiente estequiométrico do metabólito na reação atual
                    coeficiente = reacao.metabolites.get(metabolito, 0)
                    if coeficiente != 0:  # Armazena apenas valores não nulos
                    	data.append(coeficiente) # Adiciona o coeficiente à lista de dados
                    	row_ind.append(i)
                    	col_ind.append(j)

            # Criação da matriz estequiométrica esparsa no formato Compressed Sparse Row (CSR)
            matriz_S = csr_matrix((data, (row_ind, col_ind)), shape=(len(metabolitos), len(reacoes)))
            fileLog.write('MATRIZ ESTEQUIOMÉTRICA (S) - TAMANHO: {0} x {1}\n'.format(len(metabolitos), len(reacoes)))
            fileLog.write('MATRIZ S É ESPARSA: {0} elementos não nulos\n'.format(matriz_S.nnz))


		# 2. Definição da Função Objetivo
		# A função objetivo é a expressão a ser maximizada na resolução do problema de programação linear.
		# Geralmente, a função objetivo é a reação de biomassa.

		# Verifica se a função objetivo é a reação de biomassa
            biomassa_reaction = None
            # Itera sobre todas as reações no modelo para identificar a reação de biomassa
            for rxn in model.reactions:
                    # Verifica se a palavra 'biomassa' está no ID da reação (ignorando maiúsculas/minúsculas)
                    if 'biomass' in rxn.id.lower():  # Supondo que a reação de biomassa contenha 'biomass' no ID
                        biomassa_reaction = rxn
                        break
                        
            # Verifica se a reação de biomassa foi identificada
            if biomassa_reaction:
                # Define a reação de biomassa como a função objetivo do modelo (maximização)
                model.objective = {biomassa_reaction: 1.0}
                # Registra no log o ID da reação de biomassa identificada
                fileLog.write('REACAO DE BIOMASSA IDENTIFICADA: {0}\n'.format(biomassa_reaction.id))
                fileLog.write('A maximização desta reação simula o crescimento exponencial do organismo.\n')

            else:
                # Se nenhuma reação de biomassa for encontrada, registra um aviso no log
                fileLog.write('AVISO: Nenhuma reação de biomassa foi explicitamente identificada!\n')

		# Registra a função objetivo atual
            fileLog.write('FUNÇÃO OBJETIVO ATUAL: {0}\n'.format(model.objective.expression))

		# 3. Estado Pseudoestacionário
		# O FBA assume que o sistema está em estado pseudoestacionário, onde S * v = 0.


		# 4. Restrições de Capacidade (Limites Inferior e Superior)
		# As restrições de capacidade são definidas pelos limites inferior e superior dos fluxos.

		# Itera sobre todas as reações para verificar e definir limites de fluxo
            for reacao in model.reactions:
            	# Verifica se os limites de fluxo (inferior ou superior) não estão definidos
            	if reacao.lower_bound is None or reacao.upper_bound is None:
            		fileLog.write(f'AVISO: Reação {reacao.id} não possui limites de fluxo definidos. Definindo limites padrão.\n')
            		# Define limites padrão para a reação
            		reacao.lower_bound = -1000  # Limite inferior padrão
            		reacao.upper_bound = 1000   # Limite superior padrão
            		
            # Obtém os limites inferiores e superiores de todas as reações
            vlb = model.reactions.list_attr('lower_bound')
            vub = model.reactions.list_attr('upper_bound')
            
            fileLog.write('LIMITES DE FLUXO:\n')
            fileLog.write('Limite Inferior: {0}\n'.format(vlb))
            fileLog.write('Limite Superior: {0}\n'.format(vub))
            
               # 5. Otimização do Fluxo
		# O problema de FBA é formulado como:
		# Maximizar c^T * v
		# Sujeito a S * v = 0
		# E limite_inferior <= v <= limite_superior
		# Onde:
		# - c é o vetor de pesos (função objetivo, geralmente a reação de biomassa)
		# - S é a matriz estequiométrica
		# - v é o vetor de fluxos
            

		# 5. Otimização do Fluxo
		# Realiza a otimização do fluxo para maximizar a função objetivo.
            fba_solution = model.optimize()
            fileLog.write('FBA SOLUTION STATUS: {0}\n'.format(fba_solution.status))

            # Verifica se a solução encontrada é ótima
            if fba_solution.status != 'optimal':
                fileLog.write('ERRO: O modelo não encontrou uma solução ótima para FBA.\n')
                raise Exception('O modelo não encontrou uma solução ótima para FBA.')

            # Obtém o vetor de fluxos ótimo
            v_star = fba_solution.fluxes
            v_star_abs = np.abs(v_star)
            
            # Separa os fluxos em diretos e reversos
            v_2m_star = np.concatenate([
                (v_star_abs + v_star) / 2,  # Fluxos diretos
                (v_star_abs - v_star) / 2   # Fluxos reversos
            ])

		# 6. Normalização dos Fluxos pelo Peso da Biomassa Celular Seca (DCW)
		# Os fluxos são normalizados pelo peso da biomassa celular seca (DCW), com unidades de mmol/g DCW/h.

            if biomassa_reaction:
                # Obtém o fluxo de biomassa
                dcw_flux = v_star[biomassa_reaction.id]
		    
                # Normaliza os fluxos pelo fluxo de biomassa
                fluxos_normalizados = {reacao_id: fluxo / dcw_flux for reacao_id, fluxo in v_star.items()}
		    
                # Registra o fluxo de biomassa normalizado
                fileLog.write('FLUXO DE BIOMASSA (DCW NORMALIZADO): {0} mmol/g DCW/h\n'.format(dcw_flux))
            else:
                fileLog.write('AVISO: Não foi possível normalizar os fluxos, pois a reação de biomassa não foi identificada.\n')
                fluxos_normalizados = v_star  # Usa os fluxos brutos se a reação de biomassa não for encontrada


            fba_results_file = directory + "/0Marco-fba_results.csv"
            fba_results_df = pd.DataFrame({
                'Reaction': model.reactions.list_attr('id'),
                'Flux': v_star,
                'Flux Normalizado (mmol/g DCW/h)': [fluxos_normalizados.get(reacao.id, 0) for reacao in model.reactions],
                'Lower Bound': vlb,
                'Upper Bound': vub
            })
            fba_results_df.to_csv(fba_results_file, index=False)

            fileLog.write('GEROU FBA RESULTS FILE: {0}\n'.format(fba_results_file))
            fileLog.write('FBA ANALYSIS CONCLUÍDA COM SUCESSO\n')
            
            # 1. Cálculo dos Fluxos de Metabólitos e Pesos das Arestas
            #------------------------------------------------------------#
            

            # Inicialização das estruturas para armazenar os fluxos de metabólitos
            production_flows = {}  # Fluxo de produção de cada metabólito por reação
            consumption_flows = {}  # Fluxo de consumo de cada metabólito por reação

            for reacao in model.reactions:
                # Itera sobre os metabólitos e seus coeficientes estequiométricos na reação
                for metabolito, coeficiente in reacao.metabolites.items():
                        # Calcula o fluxo do metabólito na reação (fluxo da reação multiplicado pelo coeficiente estequiométrico)
                        fluxo = v_star[reacao.id] * coeficiente
                        if fluxo > 0:
                        	# Se o fluxo for positivo, o metabólito é produzido na reação
                        	if metabolito.id not in production_flows: 
                        		production_flows[metabolito.id] = {} # Inicializa um dicionário para o metabólito
                        	production_flows[metabolito.id][reacao.id] = fluxo # Armazena o fluxo de produção do metabólito na reação
                        	
                        elif fluxo < 0: # Se o fluxo for negativo, o metabólito é consumido na reação
                        	if metabolito.id not in consumption_flows: # Verifica se o metabólito já está no dicionário de consumo
                        		consumption_flows[metabolito.id] = {}
                        	# Armazena o fluxo de consumo do metabólito na reação (valor absoluto)
                        	consumption_flows[metabolito.id][reacao.id] = -fluxo
    

            # Cálculo dos fluxos entre reações (Flow_{i -> j}(X_k))
            mass_flows = {} # Armazena o fluxo de metabólitos entre pares de reações

            # Itera sobre os metabólitos que têm fluxos de produção
            for metabolito in production_flows:
                # Verifica se o metabólito também tem fluxos de consumo
                if metabolito in consumption_flows:
                    # Itera sobre as reações que produzem o metabólito
                    for reacao_i, flow_i in production_flows[metabolito].items():
                    	# Calcula o consumo total do metabólito (soma de todos os fluxos de consumo)
                    	total_consumption = sum(consumption_flows[metabolito].values())
                    	# Itera sobre as reações que consomem o metabólito
                    	for reacao_j, flow_j in consumption_flows[metabolito].items():
                        	# Calcula o fluxo do metabólito da reação i para a reação j
                        	flow_i_to_j = flow_i * (flow_j / total_consumption)
                        	# Verifica se o par de reações (i, j) já está no dicionário de fluxos
                        	if (reacao_i, reacao_j) not in mass_flows:
                        		mass_flows[(reacao_i, reacao_j)] = 0
                        	mass_flows[(reacao_i, reacao_j)] += flow_i_to_j

            # Cria um grafo direcionado para representar os fluxos entre reações
            G = nx.DiGraph()
            # Itera sobre os pares de reações e seus fluxos
            for (reacao_i, reacao_j), flow in mass_flows.items():
                G.add_edge(reacao_i, reacao_j, weight=flow)
                
            # Define a posição dos nós no grafo usando o layout Kamada-Kawai
            pos = nx.kamada_kawai_layout(G)
            

            # Desenho do grafo
            plt.figure(figsize=(12, 8))
            nx.draw(G, pos, with_labels=True, node_size=500, node_color='lightblue', edge_color='gray', font_size=8, arrows=True)
            plt.title("Mass Flow Graph (MFG)")

            # Salvar o gráfico em PNG
            output_file = directory + "/mass_flow_graph.png"
            plt.savefig(output_file, format="png", dpi=300)
            plt.close()

            fileLog.write(f'GEROU GRÁFICO DE FLUXO DE MASSA (MFG): {output_file}\n')
            
            
            
            
            #------------------------------------------------------------#

            
            reaction_id_to_number = {reaction.id: idx for idx, reaction in enumerate(model.reactions)}
            
            # Configura a precisão decimal para evitar arredondamentos indesejados
            getcontext().prec = 10

            # Geração do arquivo CSV "testefinalMarco.csv" com as informações dos nós e reações
            output_csv_file = directory + "/testefinalMarco.csv"

            # Função para formatar o peso no formato desejado
            def format_weight(weight):
            	weight_decimal = Decimal(str(weight))  # Converte para Decimal para maior precisão
            	if weight_decimal < Decimal('1e-3'):  # Usa notação científica para valores pequenos
            		return f"{weight_decimal:.2E}"
            	else:  # Usa formato decimal para valores maiores
            		return f"{weight_decimal:.9f}".rstrip('0').rstrip('.')  # Remove zeros à direita

            # Abre o arquivo CSV para escrita
            with open(output_csv_file, mode='w', newline='') as file:
            	writer = csv.writer(file, delimiter='\t')
    
            	# Escreve o cabeçalho do arquivo CSV
            	writer.writerow(["Source", "Target", "Weight"])
    
            	# Itera sobre as arestas do grafo G e escreve no arquivo CSV
            	for edge in G.edges(data=True):
            		source_id = edge[0]
            		target_id = edge[1]
            		weight = edge[2]['weight']
            		formatted_weight = format_weight(weight)  # Formata o peso
            		source_number = reaction_id_to_number[source_id]
            		target_number = reaction_id_to_number[target_id]
        
            		# Escreve os números no lugar dos nomes
            		writer.writerow([source_number, target_number, formatted_weight])

            # Registra a geração do arquivo CSV no log
            fileLog.write(f'GEROU GRÁFICO DE FLUXO DE MASSA (MFG): {output_file}\n')
            fileLog.write(f'GEROU ARQUIVO CSV DO MFG: {output_csv_file}\n')
            
            #------------------------------------------------------------#
            
            # Cálculo da betweenness centrality para cada nó no grafo G
            betweenness_centrality = nx.betweenness_centrality(G, weight='weight')  # Usa os pesos das arestas


            # Criação do arquivo "iML1515_mfg_nodes_ess-label_fba_pred.csv"
            output_nodes_file = directory + "/iML1515_mfg_nodes_ess-label_fba_pred.csv"

            # Extrai todos os nós únicos presentes no arquivo "testefinalMarco.csv"
            unique_nodes = set() # Usa um conjunto para armazenar nós únicos
            with open(output_csv_file, mode='r') as file:
            	reader = csv.reader(file, delimiter='\t') # Abre o arquivo CSV com delimitador de tabulação
            	next(reader)  # Pula o cabeçalho
            	for row in reader:
            		unique_nodes.add(int(row[0]))  # Source
            		unique_nodes.add(int(row[1]))  # Target

            # Ordena os nós em ordem crescente
            unique_nodes = sorted(unique_nodes)
            
            # Calcula o PageRank para cada nó no grafo G
            pagerank = nx.pagerank(G, weight='weight')  # Usa os pesos das arestas para calcular o PageRank
            
            # Obtém os valores do PageRank para calcular os percentis
            pagerank_values = list(pagerank.values()) # Converte os valores do PageRank em uma lista
            
            # Função para calcular o percentil de um valor em relação a uma lista de valores
            def calculate_percentile(value, values):
            	# Calcula o percentil como a porcentagem de valores na lista que são menores ou iguais ao valor dado
            	return (sum(v <= value for v in values) / len(values)) * 100

            # Encontra o valor máximo absoluto dos fluxos no vetor de fluxos ótimos (v_star)
            max_flux = np.max(np.abs(v_star))

            # Função para calcular a essencialidade de um nó usando regressão
            def calculate_essentiality_regression(node_name, v_star, max_flux):
            	# Verifica se o nó está presente no vetor de fluxos ótimos
            	if node_name in v_star:
            		
            		# Normaliza o fluxo pelo máximo fluxo absoluto
            		normalized_flux = v_star[node_name] / max_flux
        
            		# A essencialidade é o valor normalizado, podendo variar entre -1.0 e 1.0
            		essentiality_regression = normalized_flux
        
            		# Garantimos que o valor esteja dentro do intervalo [-1.0, 1.0]
            		#ssentiality_regression = max(-1.0, min(1.0, essentiality_regression))
            		return essentiality_regression
            	else:
            		return 0.0	


            # Gera o arquivo CSV com os nós e seus respectivos nomes
            with open(output_nodes_file, mode='w', newline='') as file:
            	writer = csv.writer(file, delimiter='\t')
    
            	# Escreve o cabeçalho do arquivo CSV
            	writer.writerow(["", "id", "Node Name", "PageRank", "PageRank Percentile", "Betweenness Centrality", "essentiality_regression"])
            	
            	# Itera sobre os nós únicos e escreve no arquivo CSV
            	for line_number, node_id in enumerate(unique_nodes, start=1):
                	
                	# Encontra o nome do nó correspondente ao ID
                	node_name = next((key for key, value in reaction_id_to_number.items() if value == node_id), None)
                	if node_name is None:
                		node_name = f"Unknown_{node_id}"
                	# Obtém o valor de PageRank para o nó
                	node_pagerank = pagerank.get(node_name, 0.0)
                	
                	# Calcula o percentil do PageRank do nó em relação aos demais
                	pagerank_percentile = calculate_percentile(node_pagerank, pagerank_values)
                	
                	# Obtém o valor de betweenness centrality para o nó
                	node_betweenness = betweenness_centrality.get(node_name, 0.0)
                	
                	
                	# Calcula a essencialidade contínua (regressão)
                	essentiality_regression = calculate_essentiality_regression(node_name, v_star, max_flux)



                	# Escreve a linha no arquivo CSV
                	writer.writerow([line_number, node_id, node_name, node_pagerank, pagerank_percentile, node_betweenness, essentiality_regression])
                	
                	
            #Criação da coluna Fba_pred
            # 1. Realizar a FBA original para obter o fluxo de biomassa no estado ótimo
            original_fba_solution = model.optimize()
            original_biomass_flux = original_fba_solution.fluxes[biomassa_reaction.id]

            # 2. Função para determinar a essencialidade de uma reação
            def calculate_fba_pred(reaction_id):
            	# Cria uma cópia do modelo para evitar alterações no modelo original
            	model_copy = model.copy()
    
            	# Fixa o fluxo da reação em 0 (remove a reação)
            	model_copy.reactions.get_by_id(reaction_id).bounds = (0, 0)
    
            	# Realiza a FBA no modelo modificado
            	fba_solution = model_copy.optimize()
            	
            	# Verifica se a otimização foi bem-sucedida
            	if fba_solution.status != 'optimal':
            		return 1.0  # Assume essencialidade se a otimização falhar
    
            	# Obtém o fluxo de biomassa após a remoção da reação
            	biomass_flux_after_removal = fba_solution.fluxes[biomassa_reaction.id]
    
            	# Determina a essencialidade
            	if original_biomass_flux == 0:  # Considera essencial se o fluxo de biomassa for próximo de 0
            		return 0.0
            	fba_pred = (original_biomass_flux - biomass_flux_after_removal) / original_biomass_flux
            	
            	# Garante que o valor esteja no intervalo [-1, 1]
            	fba_pred = max(-1.0, min(1.0, fba_pred))
            	
            	# Classificação binária para essencialidade
            	if fba_pred >= 0.99:  # Essencial se a biomassa cair drasticamente
            		return 1.0  
            	elif fba_pred <= 0.01:  # Não essencial se a remoção não afetar
            		return -1.0  
            	else:
            		return fba_pred  # Mantém valores intermediários para análise refinada
            	
            	return fba_pred
            	

            # 3. Aplicar a função para todas as reações e armazenar os resultados
            fba_pred_results = {reaction.id: calculate_fba_pred(reaction.id) for reaction in model.reactions}

            # 4. Atualizar o arquivo CSV com a coluna "fba_pred"
            output_nodes_file = directory + "/iML1515_mfg_nodes_ess-label_fba_pred.csv"

            with open(output_nodes_file, mode='r') as file:
            	lines = file.readlines()

            # Adiciona a coluna "fba_pred" ao cabeçalho
            header = lines[0].strip() + "\tFba_predicao\n"
            new_lines = [header]

            # Adiciona os valores de "fba_pred" para cada nó
            for line in lines[1:]:
                columns = line.strip().split('\t')
                node_id = columns[1]  # Obtém o ID do nó
                reaction_id = next((key for key, value in reaction_id_to_number.items() if value == int(node_id)), None)
                if reaction_id:
                	fba_pred = fba_pred_results.get(reaction_id, 0.0)  # Obtém o valor de fba_pred
                	# Adiciona o valor de fba_pred à linha
                	new_line = '\t'.join(columns) + f"\t{fba_pred}\n"
                	new_lines.append(new_line)
                
            # Reescreve o arquivo com a nova coluna
            with open(output_nodes_file, mode='w') as file:
            	file.writelines(new_lines)

            # Agora, corrigir o arquivo para remover a coluna "fba_pred" e renomear a coluna sem título
            with open(output_nodes_file, mode='r') as file:
            	lines = file.readlines()


            # Agora, corrigir o arquivo para remover a antiga coluna "fba_pred" e manter apenas "Fba_predicao"
            with open(output_nodes_file, mode='r') as file:
            	lines = file.readlines()

            # Identifica a posição da coluna "fba_pred" e mantém apenas a última coluna renomeada
            header_columns = lines[0].strip().split('\t')
            fba_pred_index = header_columns.index("Fba_predicao")

            # Processa as linhas para garantir a formatação correta
            corrected_lines = []
            for i, line in enumerate(lines):
            	columns = line.strip().split('\t')
    
            	# Mantém apenas a última coluna correspondente
            	corrected_line = '\t'.join(columns[:fba_pred_index] + [columns[-1]]) + "\n"
            	corrected_lines.append(corrected_line)

            # Reescreve o arquivo com a correção
            with open(output_nodes_file, mode='w') as file:
            	file.writelines(corrected_lines)

            fileLog.write(f'ATUALIZOU ARQUIVO CSV DE NÓS COM Fba_predicao: {output_nodes_file}\n')



            #####################################################
            # ITEM 3 - SIMULATION OF SINGLE KNOCKOUT OF REACTION
            #####################################################
            fileLog.write("INICIO DO ITEM 3\n")
            contador = 0
            lbAux = ubAux = 0

            reacoesZerandoFBA = []
            
            file_compound_reaction_ko = open(directory+"//react_biomass_zero_sbml_no_genes.txt", 'w', encoding="ISO-8859-1")
            fileLog.write('len(listReactPerturb) = {0}\n'.format(len(listReactPerturb)))
            for i in listReactPerturb:
                reacao = model.reactions.get_by_id(i)
                contador += 1
                
                #print(reacao)
                #print("Entrei no FOR for i in listReactPerturb:")
               

                if str(model) == "MODEL1507180020":
                    if i == "EX_fe2_e" or i == "FE2abc":
                        #print("reacao dando problema no arquivo da PAO1 2008")
                        #print("-----------------------------------")
                        continue
                
                FBABeforeDelete = model.optimize()
                
                valFBeforeDelete = round(float(FBABeforeDelete.objective_value), numCasasDec)
                if valFBeforeDelete != valInitialSolutionFmt:
                    raise Exception("1o IF - Erro! FBA nao bate com a primeira execucao")
                #print(valFBeforeDelete)
                
                #print(reacao.lower_bound)
                #print(reacao.upper_bound)
                
                lbAux = reacao.lower_bound
                ubAux = reacao.upper_bound
                reacao.lower_bound = 0
                reacao.upper_bound = 0

                #print(reacao.lower_bound)
                #print(reacao.upper_bound)
                
                FBAAfterDelete = model.optimize()
                #print(FBAAfterDelete)

                valFAfterDelete = round(float(FBAAfterDelete.objective_value), numCasasDec)
                #print(valFAfterDelete)
                
                if valFAfterDelete == 0.0 or valFAfterDelete == -0.0:
                    #print("vai gravar no arquivo")
                    #print("{0}\n".format(reacao.build_reaction_string(reacao.reaction)))
                    #file_compound_reaction_ko.write("{0}\n".format(reacao.build_reaction_string(reacao.reaction)))
                    file_compound_reaction_ko.write("{0}\n".format(reacao.build_reaction_string()))
                    reacoesZerandoFBA.append(reacao)

                reacao.upper_bound = ubAux
                reacao.lower_bound = lbAux
                lbAux = 0
                ubAux = 0

                FBAAfterRestoreModel = model.optimize()
                valFFormatAfterRestore = round(float(FBAAfterRestoreModel.objective_value), numCasasDec)
                if valFFormatAfterRestore != valInitialSolutionFmt:
                    raise Exception("2o IF - Erro! FBA nao bate com a primeira execucao")
                
                #print("-------------------------------------------------------------------")
            
            #print("Total de reacoes cujo FBA zerou = " + str(len(reacoesZerandoFBA)))
            fileLog.write('Total de reacoes cujo FBA zerou = {0}\n'.format(str(len(reacoesZerandoFBA))))
            file_compound_reaction_ko.close()
            
            #print("fechou o arquivo")
            #print("{0}".format(len(model.genes)))
            #print("{0}".format(len(model.genes) > 0))
            
            # AQUI ACONTECE UMA BIFURCACAO NA APLICACAO, ONDE ELE FAZ ESSES PASSOS CASO TENHAMOS GENES PARA
            # TRABALHAR, CASO CONTRARIO, ELE PULA ISSO TUDO, FAZ O CENARIO DESCRITO NO ELSE E SEGUE O FLUXO NO PASSO 7
            fileLog.write("INICIO DA BIFURCACAO\n")
            fileLog.write("if len(model.genes) = {0}\n".format(len(model.genes) > 0))
            if len(model.genes) > 0:
                
                #print("Extracao dos genes referentes as reacoes cujo FBA zerou - inicio")
                reacoesZerandoFBAArq = open(directory+"//"+"03-relacao_nocaute_reacao-gene.txt", "w")

                # Aqui pegamos as reacoes e extraimos os genes envolvidos com cada uma
                listaGenesAlvos = []
                for reacao in reacoesZerandoFBA:
                    reacoesZerandoFBAArq.write(str(reacao) + '\n')
                    for gene in reacao.genes:
                        #print(gene)
                        listaGenesAlvos.append(str(gene))
                        reacoesZerandoFBAArq.write(str(gene) + '\n')
                    #print("-----------------")
                    reacoesZerandoFBAArq.write('----------------\n\n')
                reacoesZerandoFBAArq.close()
                
                #print('fechou arquivo reacoesZerandoFBAArq')
                
                # Aqui com a lista de genes gerada a partir das reacoes, ordenamos e colocamos em um arquivo
                listaGenesAlvosSet = list(set(listaGenesAlvos))
                listaGenesAlvosSet.sort()
                arqGenesAlvos = open(directory+"//"+"04-genesAlvos.txt", "w")
                contador = 0
                for gene in listaGenesAlvosSet:
                    #print(gene)
                    contador += 1
                    arqGenesAlvos.write(str(contador) + " - " + str(gene) + "\n")
                arqGenesAlvos.close()
                #print("Total de genes encontrados por reacao = " + str(len(listaGenesAlvosSet)))
                #print("Extracao dos genes referentes as reacoes cujo FBA zerou - final")

                #print("GENERATED FILE %s" % reacoesZerandoFBAArq.name)
                #print("GENERATED FILE %s" % arqGenesAlvos.name)
                
               
              

                ##################################################
                # ITEM 4 - SIMULATION OF SINGLE KNOCKOUT OF GENES
                ##################################################
                fo = open(directory+"//"+"05-relacao_nocaute_gene-reacao.txt", "w")
                contador = 0
                listaGenesAlvosSet02 = []
                for i in model.genes:
                    FBABeforeDelete = model.optimize()
                    valFBeforeDelete = round(float(FBABeforeDelete.objective_value), numCasasDec)
                    if valFBeforeDelete != valInitialSolutionFmt:
                        raise Exception("1o IF - Erro! FBA nao bate com a primeira execucao")

                    cobra.manipulation.delete_model_genes(model, [str(i)])

                    FBAAfterDelete = model.optimize()
                    valFAfterDelete = round(float(FBAAfterDelete.objective_value), numCasasDec)
                    if valFAfterDelete == 0.0:
                        contador += 1
                        # reacoes_assoc = (l.reaction for l in i.reactions)
                        fo.write(str(contador) + ". Gene " + str(i) + " inibido.\n")
                        for x in i.reactions:
                            fo.write("%s : %s" % (x.id, x.reaction))
                        # fo.write("Reacoes associadas:\n%s" % ("{\n" + ",\n".join(reacoes_assoc) + "\n}"))
                        fo.write('\n\n')
                        listaGenesAlvosSet02.append(str(i))

                    cobra.manipulation.undelete_model_genes(model)

                    FBAAfterRestoreModel = model.optimize()
                    valFFormatAfterRestore = round(float(FBAAfterRestoreModel.objective_value), numCasasDec)
                    if valFFormatAfterRestore != valInitialSolutionFmt:
                        raise Exception("2o IF -Erro! FBA nao bate com a primeira execucao")

                fo.write("Total de genes inibidos que interrompem a geracao de biomassa = " + str(contador))
                fo.close()

                #print("Total de genes encontrados por nocaute direto = " + str(len(listaGenesAlvosSet02)))
                #print("GENERATED FILE %s" % fo.name)


                ############################################################
                # ITEM 5 - VERIFICATION AND UNIFICATION OF KNOCKOUT RESULTS
                ############################################################
                count = 0
                deParaGenes = open(directory+"//"+"06-deparaGenes.txt", "w")
                for geneAlvo in listaGenesAlvosSet:
                    if geneAlvo in listaGenesAlvosSet02:
                        count += 1
                        deParaGenes.write(geneAlvo)
                        deParaGenes.write('\n')
                deParaGenes.close()
                #print("Total de genes encontrados apos de/para = " + str(count))
                #print("GENERATED FILE %s" % deParaGenes.name)


                ##############################################################
                # ITEM 6 - SEARCH FOR THE CORRESPONDING EC NUMBER OF THE GENE
                ##############################################################
                k = KEGG(cache=True)
                k.TIMEOUT = TIMEOUT_SECONDS

                contador = 0
                with open(directory+"//"+"06-deparaGenes.txt", "r") as infile:
                    data = infile.read()
                my_list_genes = data.splitlines()

                fileGenesWithEC = open(directory+"//"+"07-genes_ECNumbers.txt", "w")
                fileAssocGenesEC = open(directory+"//"+"07-1-assoc_genes_ECNumbers.txt", "w")

                # AQUI ELE BUSCA TODOS OS ACRONIMOS DO KEGG RELACIONADOS COM O ORGANISMO SELECIONADO NA TELA
                df_list_organism = pd.read_excel(dir_path+"//list_organism_kegg.xls")

                df_list_organism_filter = df_list_organism[df_list_organism['name_organism'].str.contains(organismParam)]
                #print("df_list_organism_filter:  ", df_list_organism_filter)
                list_acron_kegg = df_list_organism_filter['acron_organism'].tolist()
                #print("list_acron_kegg:  ", list_acron_kegg)
                
                for acron in list_acron_kegg:
                    for gene in my_list_genes:
                        '''source_original = acron + ":" + gene
                        gene_with_underscore = re.sub(r'(?<=\D)(?=\d)|(?<=\d)(?=\D)', '_', gene)
                        source_with_underscore = acron + ":" + gene_with_underscore
                        sources = [source_original, source_with_underscore]'''
                        gene_with_underscore = re.sub(r'(?<=\D)(?=\d)|(?<=\d)(?=\D)', '_', gene)
                        sources = [
                                acron + ":" + gene,
                                acron + ":" + gene_with_underscore]
                        
                        for source in sources:
                            print("source:  ", source)
                            ec = k.link("enzyme", source)
                            print("ec  :", ec)
                                                  
                            returnKeggEC = str(ec).split() # split por espaco para separar o que veio do unicode
                            if (len(returnKeggEC)) == 0 :
                                continue
                            contador += 1
                            
                            gene_in_kegg = source.split(":")[1]
                            for ecnumber in returnKeggEC:
                                if "ec:" in ecnumber: # caso tenha o inicio "ec:"
                                	ecnumbersplitFinal = ecnumber.split(":") # novo split para prevalecer apenas o numero
                                	fileGenesWithEC.write(ecnumbersplitFinal[-1])
                                	fileGenesWithEC.write('\n')
                                	fileAssocGenesEC.write('{0};{1}\n'.format(gene_in_kegg, ecnumbersplitFinal[-1]))

                fileGenesWithEC.close()
                fileAssocGenesEC.close()

                #print("GENERATED FILE %s" % fileGenesWithEC.name)
                #print("GENERATED FILE %s" % fileAssocGenesEC.name)

            else:
                self.alternativeStepToGetECNumberWithoutGenes(directory)
            
            
            #######################################################
            # ITEM 7 - SEARCH FOR PROTEIN BY EC NUMBER NO DRUGBANK
            #######################################################
            fileLog.write("INICIO DO ITEM 7\n")
            with open(directory+"//"+"07-genes_ECNumbers.txt", "r") as infile:
                data = infile.read()
            my_list_ecnumbers = data.splitlines()
            my_list_ecnumbers = list(set(my_list_ecnumbers))
            my_list_ecnumbers.sort()


            # my_list_ecnumbers = ['2.5.1.7']

            # Aqui filtra os ECs encontrados no drugbank e simula a navegacao na pagina
            fileFilterEC = open(directory+"//"+"08-filter_ECNumbers_drugbank.txt", "w")

            for ec in my_list_ecnumbers:
                #print(ec)
                time.sleep(1) # esperar 1 segundo para acessar o link

                link01 = "https://go.drugbank.com/unearth/q?searcher=bio_entities&query="+str(ec)+"&approved=1&nutraceutical=1&illicit=1&investigational=1&withdrawn=1&experimental=1&us=0&ca=0&eu=0&commit=Apply+Filter"
                r = requests.get(link01)

                # verificar se houve retorno da pagina acessada
                if r.status_code == 200:

                    # Instancia de bs4 do resultado da 1a tela do drugbank
                    soup = BeautifulSoup(r.text, "lxml")

                    # filtro para a lista de todos os itens encontrados na busca do requests
                    list_found_ec = soup.findAll(attrs={"class": "search-result p-2 mb-4 p-sm-3 mb-sm-2"})

                    # caso tenha item, inicia as buscas para pegar o resto da informacao
                    if len(list_found_ec) > 0:

                        # itera na lista de itens encontrados no drugbank
                        for item in list_found_ec:
                            em = item.find('em')
                            ec_encontrado = em.text

                            #print(ec_encontrado, ec, ec_encontrado == ec)
                            if ec_encontrado != ec:
                                continue

                            # busca pelo link para acesso a 2a tela do drugbank
                            var = item.find('a')
                            # Montagem do link e acesso a 2a tela do drugbank para buscar as informacoes necessarias
                            link02 = "https://go.drugbank.com/"+str(var['href'])

                     
                            #r2 = requests.get(link02)
                            r2 = requests.get(link02)
                            html_content = r2.text
                            if r2.status_code == 200:

                                # Nova instancia de bs4 com as informacoes da 2a tela
                                # e pegando nome da proteina, organismo, uniprot id e drugbank id
                                soup2 = BeautifulSoup(html_content,'html.parser')
                                links = soup2.find_all('a', href=re.compile(r'http://www.uniprot.org/uniprot/'))

	                        #MarcoMexeuAqui.
                                uniprot_id = ""
                                for link_element in links:
                                	uniprot_id = link_element.text
                               
                	
                                organism_data_class = soup2.find(attrs={"class": "content-container"})
                                protein_name = organism_data_class.findAll('dd')[0].text
                                organism_name = organism_data_class.findAll('dd')[2].text
                                drugbank_id = str(var['href']).split("/")[-1]

                                fileFilterEC.write('{0};{1};{2};{3};{4}\n'.format(ec, protein_name, organism_name, uniprot_id, drugbank_id))

                            else:
                                fileFilterEC.write('{0};{1}\n'.format(ec, r.raise_for_status()))

                else:
                    fileFilterEC.write('{0};{1}\n'.format(ec, r.raise_for_status()))

            fileFilterEC.close()

            #print("GENERATED FILE %s" % fileFilterEC.name)
            
            ################################################
            # ITEM 8 - SEARCH FOR THE HOMOLOGUES AT UNIPROT
            ################################################
            fileLog.write("INICIO DO ITEM 8\n")
            print("INICIO DO ITEM 8")
            with open(directory+"//"+"08-filter_ECNumbers_drugbank.txt", "r") as infile:
                data = infile.read()
            my_list_uniprotid = data.splitlines()
            my_list_uniprotid = list(set(my_list_uniprotid))
            my_list_uniprotid = [item for item in my_list_uniprotid if len(item.split(";")) > 3]
            my_list_uniprotid.sort()
            
            
            with open(directory + "//" + "08-filter_ECNumbers_drugbank.txt", "w") as outfile:
                for item in my_list_uniprotid:
                        outfile.write(item + "\n")

            # Instancia de chamada ao UniProt
            u = UniProt()
            s = NCBIblast()

            u.TIMEOUT = TIMEOUT_SECONDS
            s.TIMEOUT = TIMEOUT_SECONDS
          

            #my_list_uniprotid = ["P42898"]
            
            # GET ALL JOBID FROM UNIPROT FOUND
            file_jobid = open(directory+"//"+"list_jobid.txt", "w", encoding="ISO-8859-1")
            start_time_first_req = time.time()

            for uniprotid in my_list_uniprotid:

                uniprotid = uniprotid.split(';')[3]

                findFasta = u.retrieve(uniprotid, "fasta")
                sequence = ""
                if isinstance(findFasta, str):
                    if (findFasta != ""):
                        sequence = findFasta
                    else:
                        continue
                
                elif isinstance(findFasta, Response):
                    if findFasta.ok:
                        sequence = findFasta.text
                    else:
                        continue
                    
                else:
                    continue
                    
             
                sequence = sequence.split("\n", 1)[1]
                sequence = re.sub(r"\n", r"", sequence)
                
               
                # executando o job que faz o blast
                jobid = s.run(program="blastp", database="uniprotkb_reference_proteomes", sequence=sequence, stype="protein", email="thiago.merigueti@ioc.fiocruz.br", alignments='1000')
                file_jobid.write("{0};{1}\n".format(uniprotid, jobid))
                #print("{0};{1}\n".format(uniprotid, jobid))
                
                time.sleep(1)
                
                print("----------------------------------------------")
                
            elapsed_time = time.time() - start_time_first_req
            #print("TEMPO TOTAL PRIMEIRO PASSO = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
            file_jobid.close()
            
            fileResultBlast = open(directory+"//"+"11-hitsEncontradosUniprot.txt", "w")
            
            # AFTER GENERATE ALL JOBID FROM BLAST, WE ITERATE ALL ITENS AND GET ALL XML RETURN
            with open(directory+"//"+"list_jobid.txt", "r") as infile:
                data = infile.read()
            my_list_jobid = data.splitlines()
            #my_list_jobid = list(set(my_list_jobid))

            arquivo = open(directory+"//"+"list_jobid.txt", 'r')

            eliminados = []
            for item in arquivo:
                #print(item)
                item = item.strip()
                uniprotid = item.split(';')[0]
                jobid = item.split(';')[1]

                url_status_blast = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
                url_result_blast = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"
                last_url_result_blast = "/xml"

                start_time_first_req = time.time()

                condition = True
                count_error = 0
                url_1 = url_status_blast+jobid
                count_refresh = 0
                while condition:
                    response = requests.get(url_1)
                    #print(response.text)
                    if response.text != "RUNNING":
                        if response.text == "FINISHED":
                            condition = False

                        if response.text == "ERROR" or response.text == "NOT_FOUND":
                            count_error += 1
                            #print("ERRO!", count_error)
                            if count_error == 3:
                                condition = False

                    time.sleep(3)
                    if count_refresh > 60:
                        condition = False

                    count_refresh += 1

                elapsed_time = time.time() - start_time_first_req
                #print("TEMPO TOTAL PRIMEIRO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

                start_time_second_req = time.time()

                if response.text != 'FINISHED':
                    continue

                url_2 = url_result_blast+jobid+last_url_result_blast
                #print("url_2:", url_2)
                response = requests.get(url_2)
                soup = BeautifulSoup(response.text, "lxml")
                #print("soup: ", soup)

                elapsed_time_2 = time.time() - start_time_second_req

                # busca por todos os hits encontrados
                listAllHits = soup.findAll('hit')

                # itera os hits e, caso encontre pseudomonas, guarda as informacoes dele
                fileListAllHist = open(dir_blasts+"//hit_organism_found_"+uniprotid+".txt", 'w')

                percentSimilarHuman = 0.0
                for hit in listAllHits:
                    organism = hit['description'].split("=")[1].split("GN")[0].strip()
                    if "Homo sapiens" in organism:
                        percentSimilarHuman = float(hit.find('identity').string.strip())
                        eValueHuman = float(hit.find('expectation').string.strip()) 


                        if percentSimilarHuman >= 40.0 and eValueHuman <= 10**-5:
                        	#print(f"Hit descartado: {uniprotid}, Identidade: {percentSimilarHuman}%, E-value: {eValueHuman}")
                        	if eValueHuman > 0:
                                	log_evalue = -math.log10(eValueHuman)
                                	eliminados.append((uniprotid, percentSimilarHuman, eValueHuman, log_evalue))
                                	print(f"Hit descartado: {uniprotid}, Identidade: {percentSimilarHuman}%, E-value: {eValueHuman}")
                                	break
                        	else:
                                	#print(f"Erro: eValueHuman não é positivo para {uniprotid}. Valor: {eValueHuman}")
                                	continue

                if percentSimilarHuman < 40.0:
                    for hit in listAllHits:
                        organism = hit['description'].split("=")[1].split("GN")[0].strip()
                        fileListAllHist.write("{0}\n".format(organism))            
             
                        if organismParam in organism: # se o valor da combo for encontrado no hit, armazena
                                idPAO1 = hit['ac']
                                percentBlast = hit.find('identity').string.strip()
                                eValue = hit.find('expectation').string.strip()
                                #print("percentBlast :" , percentBlast)
                                #print("percentSimilarHuman:",  percentSimilarHuman)
		               
                                #MarcoMexeuAqui.
                                link = f"https://www.uniprot.org/uniprot/{idPAO1}.xml"
                                #pint("link: ", link)
                                d = requests.get(link)
                                soup = BeautifulSoup(d.text, "html.parser")
		                
                                a1 = ' '.join(soup.gene.stripped_strings) if soup.gene else 0
                                a2 = "; ".join(comment_element.text.strip() for comment_element in soup('comment', {'type': 'pathway'})) or 0
                                a3 = "".join(comment_element.text.strip() for comment_element in soup('comment', {'type': 'function'})) or 0
                                a4 = "; ".join(comment_element.text.strip() for comment_element in soup('comment', {'type': 'catalytic activity'})) or 0
		                #a5 = "; ".join("; ".join(l.get_text(strip=True) for l in s('location')) for s in soup('subcellularlocation')) if soup('subcellularlocation') else 0
                                a5 = "; ".join(comment_element.text.strip().replace('\n', '; ') for comment_element in soup('subcellularlocation')).lstrip(';') or 0
                                a6 = ";".join(db['id'] for db in soup("dbreference", {"type": "PDBsum"})) or 0
		               
		                
                                if float(percentSimilarHuman) < float(percentBlast): 
		                
		                                    
                                	fileResultBlast.write("{0}---{1}---{2}---{3}---{4}---{5}---{6}---{7}---{8}---{9}\n".format(
		                        uniprotid, idPAO1, percentBlast, eValue, a1, a2, a3, a4, a5, a6
		                    ))

                                else: # Caso tenha hit com percent de humano maior, poe tudo zerado
                                	fileResultBlast.write("{0}---{1}---{2}---{3}---{4}---{5}---{6}---{7}---{8}---{9}\n".format(
		                        uniprotid, idPAO1, percentBlast, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
		                    ))

                fileListAllHist.close()
            
            fileResultBlast.close()
            arquivo.close()
            print ("Encerrou o número 8")
            
            todos_eliminados = True

            if todos_eliminados:
                eliminados.sort(key=lambda x: x[3], reverse=True)
                df_eliminados = pd.DataFrame(eliminados, columns=["UniProtID", "PercentIdentidadeHumano", "Evalue", "log(Evalue)"])
                df_eliminados.to_excel(directory + "//alvos_eliminados.xlsx", index=False)
                print("Planilha de alvos eliminados gerada com sucesso!")
            else:
                print("Nem todos os alvos foram eliminados, a planilha não será gerada.")


            ##print("GENERATED FILE %s" % fileListAllHist.name)
            #print("GENERATED FILE %s" % fileResultBlast.name)


            #####################################################
            # ITEM 9 - SEARCH FOR THE INHIBITORS FOR THE BACTERIA
            #####################################################
            fileLog.write("INICIO DO ITEM 9\n")
            # PEGANDO A ENTRADA DO ITEM POR MEIO DOS ARQUIVOS GERADOS E CONVERTENDO A DATAFRAME PARA FACILITAR
            my_list_uniprotid_drugbankid = []
            my_list_uniprotid_hit = []
            
            #data08 = pd.read_csv(directory + "//" + "08-filter_ECNumbers_drugbank.txt", sep=";", header=None)
            data08 = pd.DataFrame()
            try:
                data08 = pd.read_csv(directory + "//" + "08-filter_ECNumbers_drugbank.txt", sep=";", header=None)
            except:
                print('Note: file 08-filter... was empty. Skipping.')

            if not data08.empty:
                data08.columns = ["EC NUMBER", "PRODUCT", "ORGANISM NAME", "UNIPROTID", "DRUGBANKID"]
                data08 = data08.sort_values("EC NUMBER")
                my_list_uniprotid_drugbankid = data08["UNIPROTID"].tolist()
                data08.to_excel(directory + "//" + "08-filter_ECNumbers_drugbank.xlsx")

	    #MarcoMexeuAqui.
            data11 = pd.DataFrame()
            try:
                data11 = pd.read_csv(directory + "//" + "11-hitsEncontradosUniprot.txt", sep="---", header=None)
                if not data11.empty:
                    data11.columns = ["UNIPROTID", "HIT_UNIPROTID", "PERCENT_IDENT_BLAST", "EVALUE", "GENE_NAME",
		                          "PATHWAY", "FUNCTION", "CATALYTIC ACTIVITY", "LOCALIZATION", "ID PDB"]
                    data11 = data11.sort_values("UNIPROTID")
                    my_list_uniprotid_hit = data11["UNIPROTID"].tolist()
                    data11.to_excel(directory + "//" + "11-hits_Uniprot.xlsx")
            except FileNotFoundError:
           	 print('Note: file 11-hits... was empty. Skipping.')
            

            # GERACAO DOS ARQUIVOS DE SAIDA DESSE ITEM
            fileDataDrugs = open(directory + "//" + "13-list_inhibitors_per_target.txt", "w")
            fileInhibitorsDrugs = open(directory + "//" + "14-list_inhibitors_approved.txt", "w")

            for uniprot_drugbank in my_list_uniprotid_drugbankid:

                if uniprot_drugbank in my_list_uniprotid_hit:
                    data11Aux = data11[data11['UNIPROTID'] == uniprot_drugbank]

                    hit_toxicity = False
                    if 'GENE NAMES' in data11Aux.columns:
                       for item in data11Aux.itertuples():
                        validate = item._asdict().get('GENE NAMES') 
                        if validate == '0.0':
                            hit_toxicity = True
                            break

                    # CASO TENHA HIT MENOR QUE HUMANO, ESTARA CHEIO DE ZEROS NO RESULTADO
                    # SENDO ASSIM, ELE NAO DEVE SER CONSIDERADO UM BOM RESULTADO E NAO DEVE SER ARMAZENADO
                    if hit_toxicity:
                        continue


                    data08Aux = data08[data08["UNIPROTID"] == uniprot_drugbank]
                    drugbankid = data08Aux.iloc[0][4]


                    linkAccessForDrugsID = "https://go.drugbank.com/bio_entities/" + str(drugbankid)
                    #linkAccessForDrugsID = "https://go.drugbank.com/bio_entities/"+str(drugbankid)
                    r = requests.get(linkAccessForDrugsID)

                    if r.status_code == 200:
                        soup = BeautifulSoup(r.text)
                        listDrugsFound = soup.find(attrs={"class": "table table-sm table-bordered datatable dt-responsive"})
                        listDrugsFound = listDrugsFound.find('tbody')
                        listDrugsFound = listDrugsFound.findAll('tr')  # lista com todos os farmacos para o drugbankid informado

                        for item in listDrugsFound:
                            drugbank_drug_id = item('td')[0].text
                            drug_name = item('td')[1].text
                            drug_group = item('td')[2].text
                            pharma_action = item('td')[3].text
                            #MarcoMexeuAqui.
                            actions = ', '.join(action.text.strip() for action in item('td')[4].find_all('strong', class_='badge badge-secondary'))

                            # 0-uniprotid;;1-drugbankid;;2-drugbank_drug_id;;3-drug_name;;
                            # 4-drug_group;;5-pharma_action;;6-actions
                            fileDataDrugs.write("{0}#{1}#{2}#{3}#{4}#{5}#{6}\n".format(
                                uniprot_drugbank, drugbankid, drugbank_drug_id, drug_name,
                                drug_group, pharma_action, actions
                            ))

	                    #MarcoMexeuAqui.
                            if pharma_action == 'yes' and 'inhibitor' in actions:
                                # 0-uniprotid;;1-drugbankid;;2-drugbank_drug_id;;3-drug_name;;
                                # 4-drug_group;;5-pharma_action;;6-actions
                                fileInhibitorsDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                                    uniprot_drugbank, drugbankid, drugbank_drug_id, drug_name,
                                    drug_group, pharma_action, actions
                                ))

                    else:
                        continue
                        '''
                        fileDataDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                            uniprot_drugbank, drugbankid, r.raise_for_status(), r.raise_for_status(),
                            r.raise_for_status(), r.raise_for_status(), r.raise_for_status()
                        ))
                        '''

            fileDataDrugs.close()
            fileInhibitorsDrugs.close()

            data13 = pd.DataFrame()
            try:
                data13 = pd.read_csv(directory + "//" + "13-list_inhibitors_per_target.txt", sep="#", header=None)
            except:
                print('Note: file 13-list_inhi... was empty. Skipping.')

            #data13 = pd.read_csv(directory + "//" + "13-list_inhibitors_per_target.txt", sep=";", header=None)
            if not data13.empty:
                data13.columns = ["UNIPROTID", "DRUGBANKID", "DRUGBANKDRUGID", "DRUGNAME", "DRUGGROUP", "PHARMAACTION", "ACTIONS"]
                data13 = data13.sort_values("UNIPROTID")
                data13.to_excel(directory + "//" + "13-list_inhibitors_per_target.xlsx")

            data14 = pd.DataFrame()
            try:
                data14 = pd.read_csv(directory + "//" + "14-list_inhibitors_approved.txt", sep=";", header=None)
            except:
                print('Note: file 14-list_inhi... was empty. Skipping.')

            #data14 = pd.read_csv(directory + "//" + "14-list_inhibitors_approved.txt", sep=";", header=None)
            if not data14.empty:
                data14.columns = ["UNIPROTID", "DRUGBANKID", "DRUGBANKDRUGID", "DRUGNAME", "DRUGGROUP", "PHARMAACTION", "ACTIONS"]
                data14 = data14.sort_values("UNIPROTID")
                data14.to_excel(directory + "//" + "14-list_inhibitors_approved.xlsx")


            try:
                df_08 = pd.read_excel(directory + "//" + "08-filter_ECNumbers_drugbank.xlsx")
                df_11 = pd.read_excel(directory + "//" + "11-hits_Uniprot.xlsx")
                df_13 = pd.read_excel(directory + "//" + "13-list_inhibitors_per_target.xlsx")
    
    	        #MarcoMexeuAqui.
                df_merged = pd.merge(pd.merge(data08, data11, on='UNIPROTID', how='inner'), data13, on=['UNIPROTID','DRUGBANKID'], how='inner')
                df_merged = df_merged.drop_duplicates()
                df_merged = df_merged.drop(columns=['LOCALIZATION'])
                df_merged = df_merged.reset_index(drop=True)
                df_merged.to_excel(os.path.join(directory, "summary_results.xlsx"), index=False)
                
            except:
                print("UM DOS ARQUIVOS ESTA VAZIO. ARQUIVO NAO GERADO.")
                
            try:
                df_08 = pd.read_excel(directory + "//" + "08-filter_ECNumbers_drugbank.xlsx")
                df_11 = pd.read_excel(directory + "//" + "11-hits_Uniprot.xlsx")
    
                df_merged = pd.merge(df_08, df_11, on='UNIPROTID', how='inner')
                df_merged = pd.merge(df_merged, df_eliminados, left_on='UNIPROTID', right_on='UniProtID', how='left')
                df_merged = df_merged.drop_duplicates().reset_index(drop=True)
                df_merged.to_excel(directory + "//" + "merged_results.xlsx", index=False)
                
            except:
                print("UM DOS ARQUIVOS ESTA VAZIO. ARQUIVO NAO GERADO.")

            ##################################################################################
            # ITEM 10 - LAST ITENS FOR THE METHOD WHERE ZIP ALL DATAS GENERATED AND SEND MAIL
            ##################################################################################
            fileLog.write("INICIO DO ITEM 10\n")
            zip_file_report = self.zipReportsToSendMail(dir_path, dataHoraAtualFmt, directory)
            self.sendMailWithReportAttached(name, email, zip_file_report, model, method)
            fileLog.write("ACABOU!!! SEJA FELIZ!!!\n")

        except Exception as e:
            self.sendMailWithError(name, email, str(traceback.format_exception(None, e, e.__traceback__)))

        end = time.time()
        time_second = end - start
        time_minutes = time_second / 60
        time_hour = time_minutes / 60
        
        #print("TEMPO TOTAL DE EXECUCAO = {0} segundos".format(time_second))
        #print("TEMPO TOTAL DE EXECUCAO = {0} minutos".format(time_minutes))
        #print("TEMPO TOTAL DE EXECUCAO = {0} horas".format(time_hour))
        fileLog.close()

    # GENERATE ZIP FILE TO SEND MAIL TO USER
    def zipReportsToSendMail(self, dir_path, dataHoraAtualFmt, directory):
        zf = ZipFile(dir_path+"/results/"+dataHoraAtualFmt+'_results.zip', "w")
        for root, subdirs, files in os.walk(directory):
            for filename in files:
                file_extension = filename.split(".")[-1]
                if file_extension == 'xlsx':
                    zf.write(os.path.join(root, filename), filename, ZIP_DEFLATED)
        zf.close()

        return zf

    # SEND MAIL WITH RESULTS
    def sendMailWithReportAttached(self, name, email, zip_file_report, nameFile, method):
        methodSelected = ""
        if method == "1":
            methodSelected = "FBA+FVA"
        else:
            methodSelected = "Only FBA"
        
        if nameFile == "":
            nameFile = "NAME NOT FOUND"
        
        mensagem = "FINAL REPORT - FIND TARGETS WEB \n\n\n" \
                   "" \
                   "FILE: " + str(nameFile) + " | METHOD: " + str(methodSelected) +"\n" \
                   "" \
                   "YOUR ANALYSIS HAS FINISHED SUCESSFULLY. THE FOLLOWING FILES HAVE BEEN GENERATED:\n" \
                   "" \
                   "08-filter_ECNumbers_drugbank.xlsx\n" \
                   "FIELDS: EC NUMBER, PRODUCT, ORGANISM NAME, UNIPROTID, DRUGBANKID\n\n" \
                   "" \
                   "11-hits_Uniprot.xlsx\n" \
                   "FIELDS: UNIPROTID, HIT_UNIPROTID, PERCENT_IDENT_BLAST, EVALUE, GENE_NAME, PATHWAY, FUNCTION, CATALYTIC ACTIVITY, LOCALIZATION, ID PDB\n\n" \
                   "" \
                   "13-list_inhibitors_per_target.xlsx AND 14-list_inhibitors_approved.xlsx\n" \
                   "FIELDS: UNIPROTID	DRUGBANKID	DRUGBANKDRUGID	DRUGNAME	DRUGGROUP	PHARMAACTION	ACTIONS\n\n" \
                   "" \
                   "model_data.xlsx\n" \
                   "ALL GENES, REACTIONS AND METABOLITES IN THE SBML FILE USED IN ANALYSIS\n\n" \
                   "" \
                   "summary_results.xlsx\n" \
                   "FILE WITH A SUMMARY OF DATA FOUND IN 08, 11 AND 13 XLS FILES.\n\n" \

        # remetente    = 'findtargetsweb_fiocruz@hotmail.com'
        remetente = 'findtargetweb@gmail.com'
        # senha        = 'proccfiocruz1234'
        # senha = 'whzctyqjvqsxojcz'
        senha = 'ntucaeevfacubqyj'

        # Informacoes da mensagem
        destinatario = email
        assunto      = 'REPORT THERAPEUTICS TARGETS FROM YOUR NETWORK MODEL'
        
        msg = MIMEMultipart()
 
        msg['From'] = remetente
        msg['To'] = destinatario
        msg['Subject'] = assunto
         
        # Preparando a mensagem
        '''
        mensagem = '\r\n'.join([
          'From: %s' % remetente,
          'To: %s' % destinatario,
          'Subject: %s' % assunto,
          '',
          '%s' % mensagem
          ])
        '''
        mensagem = '\r\n'.join([
            '%s' % mensagem
        ])

        mensagem = mensagem.encode("UTF-8")
        
        msg.attach(MIMEText(mensagem.decode("UTF-8"), 'plain'))
        
        filename = zip_file_report.filename
        attachment = open(filename, "rb")
         
        part = MIMEBase('application', 'zip')
        part.set_payload((attachment).read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', "attachment", filename=os.path.basename(filename))
        msg.attach(part)
         
        # Enviando o email (USANDO O SMTP DO HOTMAIL PARA ENVIAR)
        # server = smtplib.SMTP("smtp.live.com: 587")
        server = smtplib.SMTP("smtp.gmail.com: 587")
        server.starttls()
        server.login(remetente,senha)
        text = msg.as_string()
        server.sendmail(remetente, destinatario, text)
        server.quit()


    # SEND MAIL WITH ERRORS!
    def sendMailWithError(self, name, email, dsc_exception):
        mensagem = "ERROR! PLEASE CONTACT ADMINISTRATOR OR TRY AGAIN\n\n"
        mensagem = mensagem + dsc_exception

        # remetente    = 'findtargetsweb_fiocruz@hotmail.com'
        remetente    = 'findtargetweb@gmail.com'
        # senha        = 'proccfiocruz1234'
        # senha        = 'whzctyqjvqsxojcz'
        senha = 'ntucaeevfacubqyj'

        # Informacoes da mensagem
        destinatario = email
        assunto = 'REPORT THERAPEUTICS TARGETS FROM YOUR NETWORK MODEL - ERROR!'

        msg = MIMEMultipart()

        msg['From'] = remetente
        msg['To'] = destinatario
        msg['Subject'] = assunto

        # Preparando a mensagem
        mensagem = '\r\n'.join([
            '%s' % mensagem
        ])

        mensagem = mensagem.encode("UTF-8")

        msg.attach(MIMEText(mensagem.decode("UTF-8"), 'plain'))

        # Enviando o email (USANDO O SMTP DO GMAIL PARA ENVIAR)
        #server = smtplib.SMTP("smtp.live.com: 587")
        server = smtplib.SMTP("smtp.gmail.com: 587")
        server.starttls()
        server.login(remetente, senha)
        text = msg.as_string()
        server.sendmail(remetente, destinatario, text)
        server.quit()


    # METODO QUE ARMAZENA EM PLANILHA OS DADOS DE GENE E REACAO DE UM MODELO AVALIADO
    # @PARAMS => model = modelo | nomeArquivoModelo = nome do arquivo sbml
    def reportModel(self, model, nomeArquivoModelo):
        workbook_dados_modelo = xlwt.Workbook()
        sheet_genes = workbook_dados_modelo.add_sheet("GENES")
        sheet_genes.write(0, 0, "ID GENE")
        sheet_genes.write(0, 1, "GENE")
        sheet_genes.write(0, 2, "REACTIONS")
        contador_genes = 1
        for gene in model.genes:
            assoc = (l.id for l in gene.reactions)
            sheet_genes.write(contador_genes, 0, gene.id)
            sheet_genes.write(contador_genes, 1, gene.name)
            sheet_genes.write(contador_genes, 2, ", ".join(assoc))
            contador_genes += 1
        
        sheet_react = workbook_dados_modelo.add_sheet("REACTIONS")
        sheet_react.write(0, 0, "ID REACTION")
        sheet_react.write(0, 1, "NAME")
        sheet_react.write(0, 2, "COMPOSITION")
        sheet_react.write(0, 3, "LIMITS INF/SUP")
        sheet_react.write(0, 4, "SUBSYSTEM")
        sheet_react.write(0, 5, "GENES")
        sheet_react.write(0, 6, "METABOLITES")
        contador_react = 1
        
        for reacao in model.reactions:
            assoc_genes = (l.id for l in reacao.genes)
            assoc_metab = (l.id for l in reacao.metabolites)
            sheet_react.write(contador_react, 0, reacao.id)
            sheet_react.write(contador_react, 1, reacao.name)
            sheet_react.write(contador_react, 2, reacao.build_reaction_string(reacao.reaction))
            sheet_react.write(contador_react, 3, str(reacao.bounds))
            sheet_react.write(contador_react, 4, reacao.subsystem)
            sheet_react.write(contador_react, 5, ", ".join(assoc_genes))
            sheet_react.write(contador_react, 6, ", ".join(assoc_metab))
            contador_react += 1
            
        workbook_dados_modelo.save(str(nomeArquivoModelo))
        
    
    
    # METHOD TO GET EC NUMBER FROM SBML WITHOUT GENES MAPPED
    def alternativeStepToGetECNumberWithoutGenes(self, directory):
        with open(directory+"//react_biomass_zero_sbml_no_genes.txt", "r") as infile:
            data = infile.read()
        my_list_compound_sbml = data.splitlines()
        
        k = KEGG(False, True)
        k.TIMEOUT = 500000
        pd.options.display.max_colwidth = 1000
        pd.options.display.max_rows = 1000
        
        file_ecs = open(directory+"//07-genes_ECNumbers.txt", "w", encoding="ISO-8859-1")
        file_ecs_compound = open(directory + "//07-1-assoc_EC_compounds.txt", "w", encoding="ISO-8859-1")
        file_comp_not_found_from_reactant_sbml = open(directory+"//idcomp_notfound_kegg.txt", "w", encoding="ISO-8859-1")
        file_reaction_not_found_from_reactant_sbml = open(directory+"//idreaction_not_found_kegg.txt", "w", encoding="ISO-8859-1")
        
        #print("PARTE 1 METODO ALTERNATIVO")
        for compound in my_list_compound_sbml:
            #print("parte1", compound)
            compound_no_stoich = re.sub("\d+\.\d+", "", compound) # retirada de todos os valores estequiometricos
            
            if " - reduced " in compound:
                compound_no_stoich = compound_no_stoich.replace(" - reduced ", " ")
            
            param_splt = "" # preparo do split para separar a composicao quimica
            if "<=>" in compound_no_stoich:
                param_splt = "<=>"
            else:
                param_splt = "-->"
            
            compound_splt = compound_no_stoich.split(param_splt) # separacao da composicao quimica pelo reactante(0) e o produto(1) 
            
            # Caso seja uma reacao de insercao ou uma reacao de excrecao, ele ignora
            if compound_splt[0].strip() == "" or compound_splt[1].strip() == "":
                continue  
            
            # TRATAMENTO DOS DADOS DO REAGENTE DA COMPOSICAO ENCONTRADA
            reactant_sbml = compound_splt[0].strip()
        
            list_id_cpd_kegg = []
            len_reactant_sbml = 0
            
            # verifica os reagentes que possuem mais de um composto envolvido
            if " + " in reactant_sbml: #(ex.: acetyl-coa + atp + bicarbonate)
                reactant_sbml_splt = reactant_sbml.split(" + ")
                len_reactant_sbml = len(reactant_sbml_splt)
                
                # iteracao dentro dos compostos dos reactants (acetyl-coa + atp + bicarbonate)
                for cpd_reactant_sbml in reactant_sbml_splt:
                    #http://rest.kegg.jp/find/compound/acetyl-coa
                    result_id_cpd = k.find("compound", cpd_reactant_sbml)
                    result_id_cpd_splt = result_id_cpd.split('\n')
                    result_id_cpd_splt = filter(None, result_id_cpd_splt)
                    result_id_cpd_splt_filter = list(result_id_cpd_splt)
                    
                    # iteracao dentro do resultado encontrado para um dos compostos do reagente 
                    # (cpd:C00024 Acetyl-CoA; Acetyl coenzyme A)
                    if len(result_id_cpd_splt_filter) > 0:
                        for result_cpd in result_id_cpd_splt_filter:
                            local_result_cpd_splt = result_cpd.split('\t')
                            if ";" in local_result_cpd_splt[1]:
                                dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                                for cpd_splt in dsc_cpd_splt:
                                    if cpd_splt.strip().lower() == cpd_reactant_sbml.strip().lower():
                                        list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                        break
                            else:
                                if cpd_reactant_sbml.lower() == local_result_cpd_splt[1].lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
            else:
                len_reactant_sbml = 1
                result_id_cpd = k.find("compound", reactant_sbml)
                result_id_cpd_splt = result_id_cpd.split('\n')
                result_id_cpd_splt = filter(None, result_id_cpd_splt)
                result_id_cpd_splt_filter = list(result_id_cpd_splt)
                
                if len(result_id_cpd_splt_filter) > 0:
                    for result_cpd in result_id_cpd_splt_filter:
                        local_result_cpd_splt = result_cpd.split('\t')
                        if ";" in local_result_cpd_splt[1]:
                            dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                            for cpd_splt in dsc_cpd_splt:
                                if cpd_splt.strip().lower() == reactant_sbml.strip().lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
                        else:
                            if reactant_sbml.strip().lower() == local_result_cpd_splt[1].strip().lower():
                                list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                break
        
            # caso nao tenha encontrado todos os compostos, ele nao segue o fluxo
            if len(list_id_cpd_kegg) != len_reactant_sbml:
                file_comp_not_found_from_reactant_sbml.write("{0}\n".format(compound))
                continue
            
            # PREPARO DA APLICACAO PARA USAR OS IDS DE COMPOSTOS E BUSCAR PELAS REACOES DE CADA UM
            reactant_sbml_in_cpd = "+".join(list_id_cpd_kegg)
            result_link_reactions_cpd = k.link("reaction", str(reactant_sbml_in_cpd))
            result_link_reactions_cpd_splt = result_link_reactions_cpd.split('\n')
            result_link_reactions_cpd_splt = filter(None, result_link_reactions_cpd_splt)
            result_link_reactions_cpd_splt_filter = list(result_link_reactions_cpd_splt)
            df_link_reaction_cpd = pd.DataFrame(columns=['id_cpd', 'id_reaction'])
            
            index_link_react_cpd = 0
            if len(result_link_reactions_cpd_splt_filter) > 0:
                for item_result_link_reactions_cpd in result_link_reactions_cpd_splt_filter:
                    local_item_result_link = item_result_link_reactions_cpd.split('\t')
                    df_link_reaction_cpd.loc[index_link_react_cpd] = [local_item_result_link[0], local_item_result_link[1]]
                    index_link_react_cpd += 1
            else:
                file_reaction_not_found_from_reactant_sbml.write("{0}\n".format(compound))
                continue
            
            df_link_reaction_cpd = df_link_reaction_cpd.groupby("id_reaction").filter(lambda x: len(x) == len(reactant_sbml_splt))
            set_id_reaction_kegg = {x.id_reaction for x in df_link_reaction_cpd.itertuples()}
            
            if len(set_id_reaction_kegg) > 0:
                for item_id_react in set_id_reaction_kegg:
                    result_ec_number = k.link("enzyme", item_id_react)
                    result_ec_number_splt = result_ec_number.split('\n')
                    result_ec_number_splt = filter(None, result_ec_number_splt)
                    result_ec_number_filter = list(result_ec_number_splt)
                    
                    for result_ec in result_ec_number_filter:
                        local_result_ec = result_ec.split('\t')
                        local_result_ec[1] = local_result_ec[1].replace("ec:", "").strip()
                        file_ecs.write("{0}\n".format(local_result_ec[1]))
            
        file_comp_not_found_from_reactant_sbml.close()
        file_reaction_not_found_from_reactant_sbml.close()
        
        # AQUI ELE INICIA A BUSCA POR PRODUTO DO QUE NAO FOI ENCONTRADO POR REAGENTE
        with open(directory+"//idcomp_notfound_kegg.txt", "r") as infile:
            data = infile.read()
        my_list_compound_sbml_not_found_from_reactant = data.splitlines()
        
        file_comp_not_found_from_product_sbml = open(directory+"//idcomp_notfound_2_kegg.txt", "w", encoding="ISO-8859-1")
        file_reaction_not_found_from_product_sbml = open(directory+"//idreaction_notfound_2_kegg.txt", "w", encoding="ISO-8859-1")
        
        #print("PARTE 2 METODO ALTERNATIVO")
        for compound in my_list_compound_sbml_not_found_from_reactant:
            #print("parte2", compound)
            compound_no_stoich = re.sub("\d+\.\d+", "", compound) # retirada de todos os valores estequiometricos
            
            if "- reduced" in compound:
                compound_no_stoich = compound_no_stoich.replace("- reduced", " ")
            
            param_splt = "" # preparo do split para separar a composicao quimica
            if "<=>" in compound_no_stoich:
                param_splt = "<=>"
            else:
                param_splt = "-->"
            
            compound_splt = compound_no_stoich.split(param_splt) # separacao da composicao quimica pelo reactante(0) e o produto(1) 
            
            # Caso seja uma reacao de insercao ou uma reacao de excrecao, ele ignora
            if compound_splt[0].strip() == "" or compound_splt[1].strip() == "":
                continue  
            
            product_sbml = compound_splt[1].strip()
        
            list_id_cpd_kegg = []
            len_product_sbml = 0
            
            # verifica os reagentes que possuem mais de um composto envolvido
            if " + " in product_sbml: #(ex.: acetyl-coa + atp + bicarbonate)
                product_sbml_splt = product_sbml.split(" + ")
                len_product_sbml = len(product_sbml_splt)
                
                # iteracao dentro dos compostos dos reactants (acetyl-coa + atp + bicarbonate)
                for cpd_product_sbml in product_sbml_splt:
                    #http://rest.kegg.jp/find/compound/acetyl-coa
                    result_id_cpd = k.find("compound", cpd_product_sbml)
                    result_id_cpd_splt = result_id_cpd.split('\n')
                    result_id_cpd_splt = filter(None, result_id_cpd_splt)
                    result_id_cpd_splt_filter = list(result_id_cpd_splt)
                    
                    # iteracao dentro do resultado encontrado para um dos compostos do reagente 
                    # (cpd:C00024 Acetyl-CoA; Acetyl coenzyme A)
                    if len(result_id_cpd_splt_filter) > 0:
                        for result_cpd in result_id_cpd_splt_filter:
                            local_result_cpd_splt = result_cpd.split('\t')
                            if ";" in local_result_cpd_splt[1]:
                                dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                                for cpd_splt in dsc_cpd_splt:
                                    if cpd_splt.strip().lower() == cpd_product_sbml.strip().lower():
                                        list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                        break
                            else:
                                if cpd_product_sbml.lower() == local_result_cpd_splt[1].lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
            else:
                len_product_sbml = 1
                result_id_cpd = k.find("compound", product_sbml)
                result_id_cpd_splt = result_id_cpd.split('\n')
                result_id_cpd_splt = filter(None, result_id_cpd_splt)
                result_id_cpd_splt_filter = list(result_id_cpd_splt)
                
                if len(result_id_cpd_splt_filter) > 0:
                    for result_cpd in result_id_cpd_splt_filter:
                        local_result_cpd_splt = result_cpd.split('\t')
                        if ";" in local_result_cpd_splt[1]:
                            dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                            for cpd_splt in dsc_cpd_splt:
                                if cpd_splt.strip().lower() == product_sbml.strip().lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
                        else:
                            if product_sbml.strip().lower() == local_result_cpd_splt[1].strip().lower():
                                list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                break
        
            # caso nao tenha encontrado todos os compostos, ele nao segue o fluxo
            if len(list_id_cpd_kegg) != len_product_sbml:
                file_comp_not_found_from_product_sbml.write("{0}\n".format(compound))
                continue
            
            # PREPARO DA APLICACAO PARA USAR OS IDS DE COMPOSTOS E BUSCAR PELAS REACOES DE CADA UM
            product_sbml_in_cpd = "+".join(list_id_cpd_kegg)
            result_link_reactions_cpd = k.link("reaction", str(product_sbml_in_cpd))
            result_link_reactions_cpd_splt = result_link_reactions_cpd.split('\n')
            result_link_reactions_cpd_splt = filter(None, result_link_reactions_cpd_splt)
            result_link_reactions_cpd_splt_filter = list(result_link_reactions_cpd_splt)
            df_link_reaction_cpd = pd.DataFrame(columns=['id_cpd', 'id_reaction'])
            
            index_link_react_cpd = 0
            if len(result_link_reactions_cpd_splt_filter) > 0:
                for item_result_link_reactions_cpd in result_link_reactions_cpd_splt_filter:
                    local_item_result_link = item_result_link_reactions_cpd.split('\t')
                    df_link_reaction_cpd.loc[index_link_react_cpd] = [local_item_result_link[0], local_item_result_link[1]]
                    index_link_react_cpd += 1
            else:
                file_reaction_not_found_from_product_sbml.write("{0}\n".format(compound))
                continue
            
            df_link_reaction_cpd = df_link_reaction_cpd.groupby("id_reaction").filter(lambda x: len(x) == len(product_sbml_splt))
            set_id_reaction_kegg = {x.id_reaction for x in df_link_reaction_cpd.itertuples()}
            
            if len(set_id_reaction_kegg) > 0:
                for item_id_prod in set_id_reaction_kegg:
                    result_ec_number = k.link("enzyme", item_id_prod)
                    result_ec_number_splt = result_ec_number.split('\n')
                    result_ec_number_splt = filter(None, result_ec_number_splt)
                    result_ec_number_filter = list(result_ec_number_splt)
                    
                    for result_ec in result_ec_number_filter:
                        local_result_ec = result_ec.split('\t')
                        local_result_ec[1] = local_result_ec[1].replace("ec:", "").strip()
                        file_ecs.write("{0}\n".format(local_result_ec[1]))
                        file_ecs_compound.write("{0};{1}\n".format(compound, local_result_ec[1]))
            
        file_ecs.close()
        file_ecs_compound.close()
        file_comp_not_found_from_product_sbml.close()
